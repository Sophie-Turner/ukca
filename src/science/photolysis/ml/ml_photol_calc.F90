! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
! Random forest which predicts photolysis rate coefficients
! On timesteps, using an input dataset passed from ml_photol_ctl.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Developer: Sophie Turner - st838@cam.ac.uk.
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_Photolysis
!
!  Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ml_photol_calc_mod

USE umPrintMgr, ONLY: umMessage, umPrint, umPrintFlush
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

! Subroutine available outside this module
PUBLIC :: ml_photol_calc

! Forest structure data from NetCDF. Shape = n trees, n nodes.
! Each array contains the relevant data for every tree in the forest.
! Left split children indices.
INTEGER(KIND=4), ALLOCATABLE, SAVE :: left(:,:) 
! Right split children indices.
INTEGER(KIND=4), ALLOCATABLE, SAVE :: right(:,:) 
! Features to check at nodes.
INTEGER(KIND=4), ALLOCATABLE, SAVE :: features(:,:) 
! Thresholds to split on at nodes.
REAL(KIND=4),    ALLOCATABLE, SAVE :: thresholds(:,:) 
! Node values for J rates. Shape = n trees, n J rates, n nodes.
REAL(KIND=4),    ALLOCATABLE, SAVE :: values(:,:,:) 
! Forest dimensions.
INTEGER(KIND=4), SAVE :: n_trees=0, n_nodes=0, n_outputs=0

CHARACTER(LEN=20), PARAMETER, PRIVATE :: ModuleName='ML_PHOTOL_CALC_MOD'

CONTAINS


SUBROUTINE ml_photol_calc(ml_inputs, n_xyz, n_rxns, n_features_in, preds)
! Photolysis routine which does the machine learning inference step.
! Called from src/science/photolysis/ml/ml_photol_ctl_mod.F90.

IMPLICIT NONE

! Data dimension sizes.
INTEGER, INTENT(IN) :: n_xyz, n_rxns, n_features_in
! Input and output ML datasets.
REAL, INTENT(IN)    :: ml_inputs(n_xyz, n_features_in)
REAL, INTENT(OUT)   :: preds(n_xyz, n_rxns)

! First time it's called?
LOGICAL, SAVE       :: first=.TRUE.
! Data holders for forest traversal.
REAL, ALLOCATABLE   :: preds_tree(:,:), preds_sum(:,:)
! Loop counters.
INTEGER             :: i
! Indices of J rates to mask to 0 at altitudes below mask_lvls.
INTEGER, ALLOCATABLE, SAVE :: mask_ids(:), mask_lvls(:)

! Dr Hook testing variables.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle 

CHARACTER(LEN=20),  PARAMETER :: RoutineName='ML_PHOTOL_CALC'

! Wrap it in 'Dr Hook' for timing tests.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This stuff only needs setting up at the beginning.
IF (first) THEN
  WRITE(umMessage, '(A)')                                              
  'I can see Sophies latest code changes.',
  CALL umPrint(umMessage,src=RoutineName)
  CALL umPrintFlush()
  ! Open the random forest NetCDF file
  ! and retrieve the structural arrays of the trees.
  CALL read_forest_file()
  ! Indices in prediction dataset of reactions to mask out at low altitudes
  ! and their cut_off altitudes for masking beneath.
  ! OCS, ISON, H2O, O2, N2O, MeCHO -> CH4, NO.
  ALLOCATE(mask_ids(7))
  ALLOCATE(mask_lvls(7))
  mask_ids = [5,8,12,19,21,24,25]
  mask_lvls = [65,62,68,63,62,70,65]
END IF

! Initialise prediction arrays.
ALLOCATE(preds_tree(n_xyz, n_rxns))
IF (.NOT. ALLOCATED(preds_sum)) ALLOCATE(preds_sum(n_xyz, n_rxns))
preds_tree = 0.0
preds_sum = 0.0
preds = 0.0

! Use each tree's data for regression for every gridbox.
DO i=1, n_trees
  CALL predict_tree(i, n_xyz, n_features_in, n_rxns,                           & 
                    mask_ids, mask_lvls, ml_inputs, preds_tree)
  ! Use this tree's decisions to build up a sum from all trees.
  preds_sum = preds_sum + preds_tree
END DO

DEALLOCATE(preds_tree)

! Forest predictions = average predictions from trees.
preds = preds_sum / n_trees

! Don't repeat unnecessary steps.
first = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ml_photol_calc


SUBROUTINE read_forest_file()

USE netcdf

IMPLICIT NONE

! NetCDF identity indices & return codes.
INTEGER(KIND=4) :: nc_id, dim_id, tree_id, nc_status
! Forest array indices.
INTEGER(KIND=4) :: var_id_left, var_id_right, var_id_features,                 &
                   var_id_thresholds, var_id_values 
! Test dimension size.
INTEGER(KIND=4) :: dim_len
! Test dimension name.
CHARACTER(LEN=10) :: dim_name
! Loop counter.
INTEGER(KIND=4) :: i
! Indices of tree groups and dimensions in NetCDF.
INTEGER(KIND=4), ALLOCATABLE :: tree_ids(:), dim_ids(:)
! NetCDF file path.
CHARACTER(LEN=*), PARAMETER :: filepath =                                      &
  '/home/users/sophie.turner.ext/random_forest_poc_days.nc'

! Dr Hook testing variables.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle 

CHARACTER(LEN=20),  PARAMETER :: RoutineName='READ_FOREST_FILE'

! Wrap it in 'Dr Hook' for timing tests.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nc_status = nf90_open(filepath, nf90_nowrite, nc_id)

! Get the number of trees in the forest. Maybe not needed?
nc_status = nf90_inq_dimid(nc_id, 'trees', dim_id)
nc_status = nf90_inquire_dimension(nc_id, dim_id, len=n_trees)
ALLOCATE(tree_ids(n_trees))

! Get the trees' indices. 
nc_status = nf90_inq_grps(nc_id, n_trees, tree_ids)

! Pick a tree to get the dims from, since all the trees are the same size.
tree_id = tree_ids(1)

! Get the trees' dimensions.
nc_status = nf90_inq_dimid(tree_id, 'nodes', dim_id)
nc_status = nf90_inquire_dimension(tree_id, dim_id, len=n_nodes)
nc_status = nf90_inq_dimid(tree_id, 'outputs', dim_id)
nc_status = nf90_inquire_dimension(tree_id, dim_id, len=n_outputs)

! Check that the dimension ordering is the right way around
! because the NetCDF was created in Python.
!ALLOCATE(dim_ids(2)) ! Each tree's values array is 2D.
!nc_status = nf90_inq_varid(tree_id, 'values', var_id_values)
!nc_status = nf90_inquire_variable(tree_id, var_id_values, dimids=dim_ids)
!DO i=1, 2
!  nc_status = nf90_inquire_dimension(tree_id, dim_ids(i),                      &
!                                     name=dim_name, len=dim_len)
!  WRITE(umMessage, '(A,I3,A,I8)')                                              &
!    'Dimension number, name and size in values:',                              &
!    i, dim_name, dim_len
!  CALL umPrint(umMessage,src=RoutineName)
!END DO
!CALL umPrintFlush()

! Allocate the trees' structural arrays. Shape n trees, n nodes.
ALLOCATE(left(n_trees, n_nodes))
ALLOCATE(right(n_trees, n_nodes))
ALLOCATE(features(n_trees, n_nodes))
ALLOCATE(thresholds(n_trees, n_nodes))
ALLOCATE(values(n_trees, n_outputs, n_nodes))

! Read data for every tree.
DO i=1, n_trees
  tree_id = tree_ids(i)

  ! Fetch the NetCDF arrays' indices.
  nc_status = nf90_inq_varid(tree_id, 'children_left', var_id_left)
  nc_status = nf90_inq_varid(tree_id, 'children_right', var_id_right)
  nc_status = nf90_inq_varid(tree_id, 'features', var_id_features)
  nc_status = nf90_inq_varid(tree_id, 'thresholds', var_id_thresholds)
  nc_status = nf90_inq_varid(tree_id, 'values', var_id_values)

  ! Read & write the NetCDF data into the arrays.
  nc_status = nf90_get_var(tree_id, var_id_left, left(i,:))
  nc_status = nf90_get_var(tree_id, var_id_right, right(i,:))
  nc_status = nf90_get_var(tree_id, var_id_features, features(i,:))
  nc_status = nf90_get_var(tree_id, var_id_thresholds, thresholds(i,:))
  nc_status = nf90_get_var(tree_id, var_id_values, values(i,:,:))

END DO

! Close the NetCDF.
nc_status = nf90_close(nc_id)

! +1 to all indices to convert from Python to Fortran style.
features = features + 1
left = left + 1
right = right + 1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_forest_file


SUBROUTINE predict_tree(i_tree, n_xyz, n_features_in, n_rxns,                  &
                        mask_ids, mask_lvls, ml_inputs, preds)
! Predict J rates using a single decision tree.

IMPLICIT NONE

! Which tree this is.
INTEGER, INTENT(IN) :: i_tree
! Data dimensions.
INTEGER, INTENT(IN) :: n_xyz, n_features_in, n_rxns
! Indices of J rates to set to 0.
INTEGER, INTENT(IN) :: mask_ids(7), mask_lvls(7)
! Input and output arrays.
REAL, INTENT(IN)    :: ml_inputs(n_xyz, n_features_in)
REAL, INTENT(OUT)   :: preds(n_xyz, n_rxns)

! Loop counters.
INTEGER :: i, i_sample
REAL :: lvl
! Keep track of where we are in the tree, and some important array values.
INTEGER :: node, feature
REAL :: threshold

! Dr Hook testing variables.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle 

CHARACTER(LEN=20),  PARAMETER :: RoutineName='PREDICT_TREE'

! Wrap it in 'Dr Hook' for timing tests.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                               &
!$OMP PRIVATE(i, i_sample, lvl, node, feature, threshold)                      &
!$OMP SHARED(n_xyz, i_tree, left, right, features, thresholds, values,         &
!$OMP        mask_ids, mask_lvls, ml_inputs, preds)
! Loop through samples.
DO i_sample=1, n_xyz
  ! Predict photolysis if there is daylight. Otherwise, skip and leave as 0.
  IF (ml_inputs(i_sample, 8) > 0.0) THEN
    ! Traverse tree.
    node = 1
    DO
      IF (left(i_tree, node) == 0) THEN
        ! Leaf node reached. Get the preds at this node.
        ! Values array shape: n trees, n rxns, n nodes.
        preds(i_sample, :) = values(i_tree, :, node)
        EXIT
      ELSE
        ! Keep traversing the tree.
        ! 2D tree array shapes: n trees, n nodes.
        feature = features(i_tree, node)
        threshold = thresholds(i_tree, node)
        IF (ml_inputs(i_sample, feature) <= threshold) THEN
          node = left(i_tree, node)
        ELSE
          node = right(i_tree, node)
        ENDIF
      ENDIF
    END DO
  END IF  
  
  ! Mask out negligible reactions below their cutoff model levels.
  lvl = ml_inputs(i_sample, 3)
  DO i=1, SIZE(mask_ids)
    IF (lvl < mask_lvls(i)) THEN
      WRITE(umMessage, '(A,I8)')                                              
      'Masking some J rates to 0 at low model levels. Prediction before and after:',
      preds(i_sample, mask_ids(i)
      CALL umPrint(umMessage,src=RoutineName)
      preds(i_sample, mask_ids(i)) = 0.0
      WRITE(umMessage, '(I8)')                                              
      preds(i_sample, mask_ids(i)
      CALL umPrint(umMessage,src=RoutineName)
      CALL umPrintFlush()
    END IF
  END DO
  
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN 
END SUBROUTINE predict_tree


END MODULE ml_photol_calc_mod

