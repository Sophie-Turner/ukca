! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
! Controller routine for predicting online photolysis rates
! using machine learning emulation of fast-jx 
! at pressures greater than 20 Pa
! and machine learning emulation of the 2D lookup table 
! at pressures less than 20 Pa.
! The emulation is done with a random forest.
!
! Part of the UKCA model, a community model supported by
! The Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Developer: Sophie Turner - st838@cam.ac.uk
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_Photolysis
!
! Code Description:
! Language:  FORTRAN 90
!
! ######################################################################
!
MODULE ml_photol_ctl_mod
IMPLICIT NONE
PRIVATE

! Subroutine available outside this module
PUBLIC :: ml_photol_ctl

CHARACTER(LEN=20), PARAMETER, PRIVATE :: ModuleName='ML_PHOTOL_CTL_MOD'

! Photolysis control routine which calls the machine learning inference step.
! Called from src/control/photolysis/interface/photol_ctl_mod.F90.
CONTAINS
SUBROUTINE ml_photol_ctl(                                                      &
  cols, rows, lvls, n_rxns_ukca,                                               &
  current_time, longitude, latitude, top,                                      &
  pressure, temperature, spec_humid, cloud,                                    &
  sw_flux_up, sw_flux_down, cos_sza,                                           &
  ukca_rates_names, j_rates) 

USE photol_config_specification_mod, ONLY: photol_config
USE level_heights_mod,       ONLY: r_theta_levels
USE conversions_mod,         ONLY: pi_over_180
USE umPrintMgr,              ONLY: umMessage, umPrint, umPrintFlush
USE ereport_mod,             ONLY: ereport
USE yomhook,                 ONLY: lhook, dr_hook
USE parkind1,                ONLY: jprb, jpim
USE ml_photol_calc_mod,      ONLY: ml_photol_calc
USE photol_fieldname_mod,    ONLY: photol_varname_len

IMPLICIT NONE

! Dimensions of UKCA domain.
! Num columns, aka row length. Corresponds to longitude. 16.
INTEGER, INTENT(IN)  :: cols
INTEGER, INTENT(IN)  :: rows   ! Num rows. Corresponds to latitude. 6.
INTEGER, INTENT(IN)  :: lvls   ! Num model levels.
! Num photolysis reactions requested for UKCA.
INTEGER, INTENT(IN)  :: n_rxns_ukca 

! Air pressure.
REAL, INTENT(IN)     :: pressure(cols,rows,lvls)
! Temperature.
REAL, INTENT(IN)     :: temperature(cols,rows,lvls)
! Specific humidity.
REAL, INTENT(IN)     :: spec_humid(:,:,:)
! Bulk cloud fraction.
REAL, INTENT(IN)     :: cloud(cols,rows,lvls)
! Longitude (degrees).
REAL, INTENT(IN)     :: longitude(cols,rows)
! Latitude (degrees).
REAL, INTENT(IN)     :: latitude(cols,rows)
! Top of model.
REAL, INTENT(IN)     :: top
! Current model time.
INTEGER, INTENT(IN)  :: current_time(7)
! Upward shortwave flux.
REAL, INTENT(IN)     :: sw_flux_up(cols,rows,lvls+1)
! Downward shortwave flux.
REAL, INTENT(IN)     :: sw_flux_down(cols,rows,lvls+1)
! Cosine of solar zenith angle from UM radiation scheme (not the one from UKCA).
REAL, INTENT(IN)     :: cos_sza(cols,rows)
! Names of photolysis rates requested for UKCA.
CHARACTER(LEN=photol_varname_len), POINTER, INTENT(IN) :: ukca_rates_names(:)
! Output photolysis rates.
REAL, INTENT(IN OUT) :: j_rates(cols,rows,lvls,n_rxns_ukca)

! Loop counters.
INTEGER                 :: lvl, row, col, pos, start, finish,                  & 
                           rxn_ukca, rxn_base
! Number of cols and rows, and day, from previous call.
INTEGER, SAVE           :: cols_last=0, rows_last=0, day_last=0
! Number of 3D lat-lon-alt points.
INTEGER, SAVE           :: n_xyz
! Number of photolysis base rates predicted by ML.
INTEGER                 :: n_rxns_base=26 
! ML input features: day, hour, lvl, lon, lat, sza, up sw flux, down sw flux, pres, temp.
INTEGER                 :: n_features_in=10
! Day number of the year.
INTEGER                 :: day
! Time in hours.
INTEGER                 :: hour
! Lat-lon dims changed size (unequal grid split domains)?
LOGICAL                 :: size_diff
! First time it's called?
LOGICAL, SAVE           :: first=.TRUE.
! Relative levels in case this is not an 85km, 85-level model run. Don't change.
REAL, ALLOCATABLE, SAVE :: rel_lvls(:)
! ML input array of samples x features. Saved because lvls don't need resetting.
REAL, ALLOCATABLE, SAVE :: ml_inputs(:,:)
! ML predictions array of samples x features.
REAL, ALLOCATABLE       :: preds(:,:)
! Base photolysis rates predicted by ML, used for other photolysis rates.
CHARACTER(LEN=6), ALLOCATABLE, SAVE :: base_rates_names(:)
! Indices for mapping the predicted J rates to the rates requested for UKCA.
INTEGER, ALLOCATABLE, SAVE          :: map_ids(:) 

! Dr Hook testing variables.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle 

CHARACTER(LEN=20),  PARAMETER :: RoutineName='ML_PHOTOL_CTL'

! Wrap it in 'Dr Hook' for timing tests.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! See if the lat-lon dimensions have changed.
size_diff = ((cols /= cols_last) .OR. (rows /= rows_last))
IF (size_diff) THEN
  ! Calculate number of 3D elements in domain.
  n_xyz = cols * rows * lvls  
  ! Adjust size of ML input array and preds.
  IF (ALLOCATED(ml_inputs)) DEALLOCATE(ml_inputs)
  ALLOCATE(ml_inputs(n_xyz, n_features_in))
  cols_last = cols
  rows_last = rows
END IF

! Preds are never reused so always need an allocation check.
IF (.NOT. ALLOCATED(preds)) ALLOCATE(preds(n_xyz, n_rxns_base))

! Levels & rates list to calc don't change during a model run.
IF (first) THEN ! First time photolysis is called in the program.  
  
  ! Levels.
  ALLOCATE(rel_lvls(lvls))
  DO lvl=1, lvls
    ! Estimate relative levels in case this isn't an 85-level, 85km model run.
    rel_lvls(lvl) = ((top / 1000.0) / REAL(lvls)) * REAL(lvl) 
  END DO
  
  ! List of J rates to return. They come out of the ML in this order.
  ALLOCATE(base_rates_names(n_rxns_base))
  base_rates_names = (/ 'jhchoa', 'jhchob', 'jmkal ', 'jcl2o2', 'jcos  ',      &
                      'jso3  ', 'jmena ', 'jiprn ', 'jaceta', 'jetcho',        &
                      'jno3a ', 'jh2o  ', 'jhobr ', 'jhocl ', 'jhono2',        &
                      'jpna33', 'jh2o2 ', 'jmhp  ', 'jo2   ', 'jo3a  ',        &
                      'jn2o  ', 'jmacr ', 'jmacro', 'jacetb', 'jno   ',        &
                      'jno2  ' /) 
  ! Map the predictions to the requested J rates and save these indices.
  ! Not all J rates are present in base rates and we want to skip those 
  ! and leave their original array values intact, so set missing indices 
  ! to 0 which means skip assignment.
  ALLOCATE(map_ids(n_rxns_ukca))
  map_ids = 0
  DO rxn_ukca=1, n_rxns_ukca
    DO rxn_base=1, n_rxns_base 
      IF (TRIM(ukca_rates_names(rxn_ukca)) ==                                  &
          TRIM(base_rates_names(rxn_base))) THEN
        map_ids(rxn_ukca) = rxn_base
        EXIT
      END IF
    END DO
  END DO
  
END IF ! First time photolysis is called in the program.

hour = current_time(4)
day = current_time(7)

! Build up the 2D array of samples x features for random forest.

! Pad time as features implicitly from ints. 
! No need to loop since this is called on timseteps.
! Keep the previous column if we're still on the same day.
IF (day /= day_last) THEN 
  ml_inputs(:,1) = day
  day_last = day
END IF  
ml_inputs(:,2) = hour

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                               &
!$OMP PRIVATE(lvl,row,col,pos,start,finish)                                    &
!$OMP SHARED(lvls,rows,cols,size_diff,rel_lvls,latitude,longitude,sw_flux_up,  &
!$OMP        sw_flux_down,pressure,temperature,cos_sza,ml_inputs)
! Transform 2D and 3D spatial-dim arrays into the 2D ML-dim array.
DO lvl=1, lvls
  ! Pad level to fill the samples in its feature column, 
  ! or keep the previous column if it doesn't need to change.
  IF (size_diff) THEN
    start = (lvl-1)*rows*cols + 1
    finish = lvl*rows*cols 
    ml_inputs(start:finish,3) = rel_lvls(lvl)
  END IF
  DO row=1, rows
    DO col=1, cols
      pos = (lvl-1)*rows*cols + (row-1)*cols + col 
      ml_inputs(pos,4) = latitude(col,row) 
      ml_inputs(pos,5) = longitude(col,row) 
      ml_inputs(pos,6) = cos_sza(col,row)     
      ml_inputs(pos,7) = sw_flux_up(col,row,lvl)
      ml_inputs(pos,8) = sw_flux_down(col,row,lvl)
      ml_inputs(pos,9) = pressure(col,row,lvl)
      ml_inputs(pos,10) = temperature(col,row,lvl)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
      
! Checks.
!WRITE(umMessage,'(A,2E12.3)') 'Relative levels (min, max):',                   &
!                 MINVAL(rel_lvls), MAXVAL(rel_lvls)
!CALL umPrint(umMessage,src=RoutineName) 
!WRITE(umMessage,'(A,2E12.3)') 'Day in ML input data (min, max):',              &
!                 MINVAL(ml_inputs(:,1)), MAXVAL(ml_inputs(:,1))
!CALL umPrint(umMessage,src=RoutineName)
!WRITE(umMessage,'(A,2E12.3)') 'Hour in ML input data (min, max):',             &
!                 MINVAL(ml_inputs(:,2)), MAXVAL(ml_inputs(:,2))
!CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,2E12.3)') 'Level in ML input data (min, max):',            &
                 MINVAL(ml_inputs(:,3)), MAXVAL(ml_inputs(:,3))
CALL umPrint(umMessage,src=RoutineName)
! 1st call lon = 0 to 28.
!WRITE(umMessage,'(A,2E12.3)') 'Lat in ML input data (min, max):',              &
!                 MINVAL(ml_inputs(:,4)), MAXVAL(ml_inputs(:,4)) 
!CALL umPrint(umMessage,src=RoutineName)
! 1st call lat = -89 to -83.
!WRITE(umMessage,'(A,2E12.3)') 'Lon in ML input data (min, max):',              &
!                 MINVAL(ml_inputs(:,5)), MAXVAL(ml_inputs(:,5))
!CALL umPrint(umMessage,src=RoutineName)
! 1st call cos(SZA) = 0.28 to 0.38 ~ 70 degrees.
!WRITE(umMessage,'(A,2E12.3)') 'SZA in ML input data (min, max):',              &
!                 MINVAL(ml_inputs(:,6)), MAXVAL(ml_inputs(:,6))
!CALL umPrint(umMessage,src=RoutineName)
! 1st call up sw flux = 0 to 1e-7.
!WRITE(umMessage,'(A,2E12.3)') 'Up sw flux in ML input data (min, max):',       &
!                 MINVAL(ml_inputs(:,7)), MAXVAL(ml_inputs(:,7))
!CALL umPrint(umMessage,src=RoutineName)
! 1st call down sw flux = 0 to 0.
!WRITE(umMessage,'(A,2E15.5)') 'Down sw flux in ML inputs (min, max):',         &
!                 MINVAL(ml_inputs(:,8)), MAXVAL(ml_inputs(:,8))
!CALL umPrint(umMessage,src=RoutineName)
!WRITE(umMessage,'(A,2E12.3)') 'Pressure in ML input data (min, max):',         &
!                 MINVAL(ml_inputs(:,9)), MAXVAL(ml_inputs(:,9))
!CALL umPrint(umMessage,src=RoutineName)
! 1st call temp = 195 to 288 K.
!WRITE(umMessage,'(A,2E12.3)') 'Temp in ML input data (min, max):',             &
!                 MINVAL(ml_inputs(:,10)), MAXVAL(ml_inputs(:,10))
!CALL umPrint(umMessage,src=RoutineName)
CALL umPrintFlush()          
      
! Pass the inputs to the random forest. Returns preds array.
! We only want base rates, not full array with copies.
CALL ml_photol_calc(ml_inputs, n_xyz, n_rxns_base, n_features_in, preds)

! Map the base rates onto the full array of J rates.
! Turn the predictions from 2D samples, features to 4D lon, lat, lvl, rxn.
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                               &
!$OMP PRIVATE(lvl,row,col,rxn_ukca,pos)                                        &
!$OMP SHARED(lvls,rows,cols,n_rxns_ukca,map_ids,j_rates,preds)
DO lvl=1, lvls
  DO row=1, rows
    DO col=1, cols
      pos = (lvl-1)*rows*cols + (row-1)*cols + col
      DO rxn_ukca=1, n_rxns_ukca
        IF (map_ids(rxn_ukca) /= 0) THEN
          j_rates(col,row,lvl,rxn_ukca) = preds(pos,map_ids(rxn_ukca)) 
        END IF
      END DO
    END DO
  END DO  
END DO
!$OMP END PARALLEL DO

! Deallocate preds in case domain changes.
IF (ALLOCATED(preds)) DEALLOCATE(preds)

! Don't repeat anything unnecessary.
first = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ml_photol_ctl
END MODULE ml_photol_ctl_mod                          
