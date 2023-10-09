#include "version.inc"

MODULE CALC_FORCE_DEM_MOD

   USE CFRELVEL_MOD, ONLY: CFRELVEL

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_FORCE_DEM                                          !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Calculate contact force and torque on particle from        !
!           particle-particle and particle-wall collisions. Treats     !
!           wall interaction also as a two-particle interaction but    !
!           accounting for the wall properties                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CALC_FORCE_DEM

! Modules
!---------------------------------------------------------------------//
      USE calc_collision_wall
      USE constant, ONLY: Pi
      USE des_thermo
      USE des_thermo_cond
      USE discretelement
      USE resize
      USE run
      use param1, only: one, small_number, zero, undefined

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0
! particle no. indices
      INTEGER :: I, LL, cc, CC_START, CC_END
! the overlap occurring between particle-particle or particle-wall
! collision in the normal direction
      DOUBLE PRECISION :: OVERLAP_N, OVERLAP_T(3)
! square root of the overlap
      DOUBLE PRECISION :: SQRT_OVERLAP
! distance vector between two particle centers or between a particle
! center and wall when the two surfaces are just at contact (i.e. no
! overlap)
      DOUBLE PRECISION :: R_LM,DIST_CI,DIST_CL
! the normal and tangential components of the translational relative
! velocity
      DOUBLE PRECISION :: V_REL_TRANS_NORM, rad
! distance vector between two particle centers or between a particle
! center and wall at current and previous time steps
      DOUBLE PRECISION :: DIST(3), NORMAL(3), DIST_MAG, POS(3)
! tangent to the plane of contact at current time step
      DOUBLE PRECISION :: VREL_T(3)
! normal and tangential forces
      DOUBLE PRECISION :: FN(3), FT(3)
! temporary storage of force
      DOUBLE PRECISION :: FC_TMP(3)
! temporary storage of force for torque
      DOUBLE PRECISION :: TOW_FORCE(3)
! temporary storage of torque
      DOUBLE PRECISION :: TOW_TMP(3,2)
! temporary storage of conduction/radiation
      DOUBLE PRECISION :: QQ_TMP

! store solids phase index of particle (i.e. pijk(np,5))
      INTEGER :: PHASEI, PHASELL
! local values used spring constants and damping coefficients
      DOUBLE PRECISION :: ETAN_DES, ETAT_DES
      DOUBLE PRECISION :: KN_DES, KT_DES
! local values used for calculating cohesive forces
      DOUBLE PRECISION :: FORCE_COH, EQ_RADIUS, DistApart

      LOGICAL, PARAMETER :: report_excess_overlap = .FALSE.

      DOUBLE PRECISION :: FNMD, FTMD, MAG_OVERLAP_T, TANGENT(3)
      integer :: start_index, end_index, do_index

! Force chain data
      INTEGER :: FCHAINC_OLD           ! old value of FCHAINC
      INTEGER :: old_size,new_size ! New force chain array size when resizing (growing) arrays
      DOUBLE PRECISION :: FC_MAG_MAX, MAG_ORIENT
      INTEGER :: FC_MAX_INDEX         ! FC index where max occurs (within LL-I pair)

! Rolling friction
      DOUBLE PRECISION :: RFTOW(3), omega_mag

! change the coefficient of restitution--yuxing
! Store the particle-pair imformation last time-step
      DOUBLE PRECISION :: DES_VEL_L_OLD(3), DES_VEL_I_OLD(3), DES_RADIUS_L_OLD, DES_RADIUS_I_OLD,DES_OMEGA_L_OLD(3), DES_OMEGA_I_OLD(3)
      DOUBLE PRECISION :: NORM_OLD(3), DIST_LI_OLD, V_REL_TRANS_NORM_OLD, VREL_T_OLD(3), DIST_OLD(3), DES_POS_L_OLD(3), DES_POS_I_OLD(3)

!......................................................................!

! Initialize cohesive forces
      IF(USE_COHESION) PostCohesive(:) = ZERO

      CALL CALC_DEM_FORCE_WITH_WALL_STL

      if (optflag1.eq.1) then
        start_index = 1
        end_index = MAX_PIP_EXIST
      else
        start_index = 1
        end_index = MAX_PIP
      endif


! Check particle LL neighbor contacts
!---------------------------------------------------------------------//

!$omp parallel default(none) private(pos,rad,cc,cc_start,cc_end,ll,i,  &
!$omp    overlap_n,vrel_t,v_rel_trans_norm,sqrt_overlap,dist,r_lm,     &
!$omp    kn_des,kt_des,phasell,phasei,etan_des,                        &
!$omp    etat_des,fn,ft,overlap_t,tangent,mag_overlap_t,               &
!$omp    eq_radius,distapart,force_coh,dist_mag,NORMAL,ftmd,fnmd,      &
!$omp    dist_cl, dist_ci, fc_tmp, tow_tmp, tow_force, qq_tmp,         &
!$omp    mew_r,rftow,omega_new,omega_mag )                             &
!!$omp    fchainc,fchainc_old,fchain_midpoint,fchain_fn,fchain_ft,      &
!!$omp    fchain_length,fchain_fn_mag,old_size,new_size,fchain_fcmax,   &
!!$omp    fchain_local_id1,fchain_local_id2,                            &
!!$omp    fchain_global_id1,fchain_global_id2,fchain_overlap )          &
!$omp shared(max_pip,neighbors,neighbor_index,des_pos_new,des_radius,  &
!$omp    des_coll_model_enum,kn,kt,pft_neighbor,pijk,hert_kn,hert_kt,  &
!$omp    des_etan,des_etat,mew,use_cohesion, calc_cond_des, dtsolid,   &
!$omp    van_der_waals,vdw_outer_cutoff,vdw_inner_cutoff,              &
!$omp    hamaker_constant,asperities,surface_energy, optflag1,         &
!$omp    start_index, end_index,                                       &
!$omp    tow, fc, energy_eq, grav_mag, postcohesive, pmass, q_source,  &
!$omp    write_force_chain, des_col_force, iglobal_id,                 &
!$omp    old_size,new_size,                                            &
!$omp    fchain_fcmax,fc_mag_max,fc_max_index, mag_orient,             &
!$omp    FCHAINC,FCHAINC_OLD,                                          &
!$omp    FCHAIN_MIDPOINT,FCHAIN_ORIENT, FCHAIN_FN,FCHAIN_FT,           &
!$omp    FCHAIN_LENGTH,FCHAIN_FN_MAG,FCHAIN_LOCAL_ID1,FCHAIN_LOCAL_ID2,&
!$omp    FCHAIN_GLOBAL_ID1,FCHAIN_GLOBAL_ID2, FCHAIN_OVERLAP) &
!$omp shared(max_pip_exist)

! Force chain initialization
      IF(WRITE_FORCE_CHAIN) THEN
         FCHAINC_OLD = min(FCHAINC, size(FCHAIN_MIDPOINT,1))
         FCHAINC = 0
         FCHAIN_MIDPOINT(:,:) = UNDEFINED
         FCHAIN_ORIENT(:,:) = UNDEFINED
         FCHAIN_FCMAX(:) = .FALSE.
      ENDIF

!$omp do reduction (+: FCHAINC)

      DO do_index = start_index, end_index
         if (optflag1.eq.1) then
            ll = do_index
         else
            ll = do_index
            IF(IS_NONEXISTENT(LL)) CYCLE
         endif

         IF(WRITE_FORCE_CHAIN) THEN
            FC_MAG_MAX = ZERO
            FC_MAX_INDEX = -1
         ENDIF


         pos = DES_POS_NEW(LL,:)
         rad = DES_RADIUS(LL)

         CC_START = 1
         IF (LL.gt.1) CC_START = NEIGHBOR_INDEX(LL-1)
         CC_END   = NEIGHBOR_INDEX(LL)

         DO CC = CC_START, CC_END-1
            I  = NEIGHBORS(CC)
            IF(IS_NONEXISTENT(I)) CYCLE

            R_LM = rad + DES_RADIUS(I)
            DIST(:) = DES_POS_NEW(I,:) - POS(:)
            DIST_MAG = dot_product(DIST,DIST)

            FC_TMP(:) = ZERO

! Compute particle-particle VDW cohesive short-range forces
            IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
               IF(DIST_MAG < (R_LM+VDW_OUTER_CUTOFF)**2) THEN
                  EQ_RADIUS = 2d0 * DES_RADIUS(LL)*DES_RADIUS(I) /     &
                    (DES_RADIUS(LL)+DES_RADIUS(I))
                  IF(DIST_MAG > (VDW_INNER_CUTOFF+R_LM)**2) THEN
                     DistApart = (SQRT(DIST_MAG)-R_LM)
                     FORCE_COH = HAMAKER_CONSTANT * EQ_RADIUS /           &
                        (12d0*DistApart**2) * (Asperities/(Asperities+    &
                        EQ_RADIUS) + ONE/(ONE+Asperities/DistApart)**2 )
                  ELSE
                     FORCE_COH = 2d0 * PI * SURFACE_ENERGY * EQ_RADIUS *  &
                       (Asperities/(Asperities+EQ_RADIUS) + ONE/          &
                       (ONE+Asperities/VDW_INNER_CUTOFF)**2 )
                  ENDIF
                  FC_TMP(:) = DIST(:)*FORCE_COH/SQRT(DIST_MAG)
                  TOW_TMP(:,:) = ZERO

! just for post-processing mag. of cohesive forces on each particle
                  PostCohesive(LL) = PostCohesive(LL) + FORCE_COH / PMASS(LL)
               ENDIF
            ENDIF

            IF (ENERGY_EQ) THEN
            ! Calculate conduction and radiation for thermodynamic neighbors
               IF(CALC_COND_DES(PIJK(LL,5))) THEN
                  QQ_TMP = DES_CONDUCTION(LL, I, sqrt(DIST_MAG), PIJK(LL,5), PIJK(LL,4))

!$omp atomic
                  Q_Source(LL) = Q_Source(LL) + QQ_TMP

!$omp atomic
                  Q_Source(I) = Q_Source(I) - QQ_TMP
               ENDIF
            ENDIF

            IF(DIST_MAG > (R_LM - SMALL_NUMBER)**2) THEN
               PFT_NEIGHBOR(:,CC) = 0.0
               CYCLE
            ENDIF

            IF(DIST_MAG == 0) THEN
               WRITE(*,8550) LL, I
               ERROR_STOP "division by zero"
 8550 FORMAT('distance between particles is zero:',2(2x,I10))
            ENDIF

            DIST_MAG = SQRT(DIST_MAG)
            NORMAL(:)= DIST(:)/DIST_MAG

! Calculate the normal overlap
            OVERLAP_N = R_LM-DIST_MAG
            IF(REPORT_EXCESS_OVERLAP) CALL PRINT_EXCESS_OVERLAP

! Calculate the components of translational relative velocity for a
! contacting particle pair and the tangent to the plane of contact
            CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, VREL_T,            &
               NORMAL(:), DIST_MAG)

            phaseLL = PIJK(LL,5)
            phaseI = PIJK(I,5)

! change the coefficient of restitution--yuxing
         ! calculate the particle-pair relative velocity last timestep: V_REL_TRANS_NORM_OLD, NORM_OLD
            DES_POS_L_OLD = DES_POS_OLD_Eu(L,:)
            DES_POS_I_OLD = DES_POS_OLD_Eu(I,:)

            DES_VEL_L_OLD = DES_VEL_OLD_Eu(L,:)
            DES_VEL_I_OLD = DES_VEL_OLD_Eu(I,:)
            DES_OMEGA_L_OLD = OMEGA_OLD_Eu(L,:)
            DES_OMEGA_I_OLD = OMEGA_OLD_Eu(I,:)

            DES_RADIUS_L_OLD = DES_RADIUS(L,:)
            DES_RADIUS_I_OLD = DES_RADIUS(I,:)
            
            DIST_OLD(:) = DES_POS_I_OLD - DES_POS_L_OLD
            DIST_LI_OLD = dot_product(DIST,DIST)
            DIST_LI_OLD = SQRT(DIST_LI_OLD)
            NORM_OLD(:)= DIST_OLD(:)/DIST_LI_OLD

            Call Cal_Rel_Vel(DES_VEL_L_OLD, DES_VEL_I_OLD, DES_RADIUS_L_OLD, DES_RADIUS_I_OLD,DES_OMEGA_L_OLD, DES_OMEGA_I_OLD, &
            NORM_OLD, DIST_LI_OLD, V_REL_TRANS_NORM_OLD, VREL_T_OLD)

         ! change the collosion paremeters based on velocity change: V_REL_TRANS_NORM_OLD, NORM_OLD, V_REL_TRANS_NORM, NORMAL(:)
            if (( V_REL_TRANS_NORM_OLD .ge. 0.0) .and. (abs(V_REL_TRANS_NORM) .gt. abs(V_REL_TRANS_NORM_OLD)))  then
               if ( V_REL_TRANS_NORM * V_REL_TRANS_NORM_OLD .ge. 0.0 ) then
                  IF(CALL_USR) then
                     CALL USR1_DES(abs(V_REL_TRANS_NORM),L,I)
                  end if
               end if
            end if

            if (( V_REL_TRANS_NORM_OLD .lt. 0.0) .and. (abs(V_REL_TRANS_NORM) .lt. abs(V_REL_TRANS_NORM_OLD)))  then
               IF(CALL_USR) then
                  CALL USR1_DES(abs(V_REL_TRANS_NORM),L,I)
               end if
            end if

! Hertz spring-dashpot contact model
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               sqrt_overlap = SQRT(OVERLAP_N)
               KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap
               KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap
               sqrt_overlap = SQRT(sqrt_overlap)
               ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap
               ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap

! Linear spring-dashpot contact model
            ELSE
               KN_DES = KN
               KT_DES = KT
               ETAN_DES = DES_ETAN(phaseLL,phaseI)
               ETAT_DES = DES_ETAT(phaseLL,phaseI)
            ENDIF

! Calculate the normal contact force
            FN(:) =  -(KN_DES * OVERLAP_N * NORMAL(:) + &
               ETAN_DES * V_REL_TRANS_NORM * NORMAL(:))

! Calculate the tangential overlap
            OVERLAP_T(:) = DTSOLID*VREL_T(:) + PFT_NEIGHBOR(:,CC)
            MAG_OVERLAP_T = sqrt(dot_product(OVERLAP_T,OVERLAP_T))

! Calculate the tangential contact force.
            IF(MAG_OVERLAP_T > 0.0) THEN
! Calculate the tangential contact force.
               FT = -KT_DES*OVERLAP_T - ETAT_DES*VREL_T
               FTMD = sqrt(dot_product(FT,FT))
! Max force before the on set of frictional slip.
               FNMD = MEW*sqrt(dot_product(FN,FN))
! Frictional slip
               IF(FTMD > FNMD) THEN
! Direction of tangential force.
                  TANGENT = OVERLAP_T/MAG_OVERLAP_T
                  FT = -FNMD * TANGENT
                  OVERLAP_T = (FNMD/KT_DES) * TANGENT
               ENDIF
            ELSE
               FT = 0.0
            ENDIF

! Save tangential displacement history
            PFT_NEIGHBOR(:,CC) = OVERLAP_T(:)

! calculate the distance from the particles' centers to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
            DIST_CL = DIST_MAG/2.d0 + (DES_RADIUS(LL)**2 - &
               DES_RADIUS(I)**2)/(2.d0*DIST_MAG)

            DIST_CI = DIST_MAG - DIST_CL

            TOW_force(:) = DES_CROSSPRDCT(NORMAL(:), FT(:))
            TOW_TMP(:,1) = DIST_CL*TOW_force(:)
            TOW_TMP(:,2) = DIST_CI*TOW_force(:)

! Calculate the total force FC of a collision pair
! total contact force ( FC_TMP may already include cohesive force)
            FC_TMP(:) = FC_TMP(:) + FN(:) + FT(:)

            FC(LL,:) = FC(LL,:) + FC_TMP(:)

!$omp atomic
            FC(I,1) = FC(I,1) - FC_TMP(1)
!$omp atomic
            FC(I,2) = FC(I,2) - FC_TMP(2)
!$omp atomic
            FC(I,3) = FC(I,3) - FC_TMP(3)

! for each particle the signs of norm and ft both flip, so add the same torque
            TOW(LL,:) = TOW(LL,:) + TOW_TMP(:,1)

!$omp atomic
            TOW(I,1)  = TOW(I,1)  + TOW_TMP(1,2)
!$omp atomic
            TOW(I,2)  = TOW(I,2)  + TOW_TMP(2,2)
!$omp atomic
            TOW(I,3)  = TOW(I,3)  + TOW_TMP(3,2)

            DES_COL_FORCE(I,:) = FC(I,:)

!######################################################################
! Force chain data - Implemented by Eric Breard (U. of Oregon)
!                    Reviewed by Jeff Dietiker
!######################################################################
! Note: Force Chain Data does not work in SMP
!######################################################################

            IF(WRITE_FORCE_CHAIN) THEN
! Resize arrays if needed
               IF(FCHAINC>=FCHAINC_OLD.AND.MOD(FCHAINC,100)==0) THEN
                  old_size = size(FCHAIN_MIDPOINT,1)
                  IF(FCHAINC>=INT(0.9*old_size)) THEN
                     new_size = INT(1.5*old_size)
                     call real_grow2_reverse(FCHAIN_MIDPOINT,new_size)
                     FCHAIN_MIDPOINT(old_size+1:new_size,:) = UNDEFINED
                     call real_grow2_reverse(FCHAIN_ORIENT,new_size)
                     call real_grow2_reverse(FCHAIN_FN,new_size)
                     call real_grow2_reverse(FCHAIN_FT,new_size)
                     call real_grow(FCHAIN_LENGTH,new_size)
                     call real_grow(FCHAIN_FN_MAG,new_size)
                     call logical_grow(FCHAIN_FCMAX,new_size)
                     call integer_grow(FCHAIN_LOCAL_ID1,new_size)
                     call integer_grow(FCHAIN_LOCAL_ID2,new_size)
                     call integer_grow(FCHAIN_GLOBAL_ID1,new_size)
                     call integer_grow(FCHAIN_GLOBAL_ID2,new_size)
                     call real_grow(FCHAIN_OVERLAP,new_size)
                  ENDIF
               ENDIF
! Save data
               FCHAINC                     = FCHAINC + 1  ! Force chain counter
               FCHAIN_MIDPOINT(FCHAINC,:)  = 0.5D0 * (DES_POS_NEW(LL,:) + DES_POS_NEW(I,:)) ! Mid point coordinates
               MAG_ORIENT                  = dsqrt(dot_product(FN,FN))
               if(MAG_ORIENT/=ZERO) FCHAIN_ORIENT(FCHAINC,:)    = FN(:)/MAG_ORIENT ! Orientation (same as Normal force)
               FCHAIN_FN(FCHAINC,:)        = FN(:)    ! Normal force
               FCHAIN_FT(FCHAINC,:)        = FT(:)    ! Tangential force
               FCHAIN_LENGTH(FCHAINC)      = DIST_MAG !distance between LL and I
               FCHAIN_FN_MAG(FCHAINC)      = DSQRT(FN(1)**2 + FN(2)**2 + FN(3)**2) ! Normal force magnitude
               FCHAIN_LOCAL_ID1(FCHAINC)   = LL
               FCHAIN_LOCAL_ID2(FCHAINC)   = I
               FCHAIN_GLOBAL_ID1(FCHAINC)  = iglobal_id(LL)
               FCHAIN_GLOBAL_ID2(FCHAINC)  = iglobal_id(I)
               FCHAIN_OVERLAP(FCHAINC)     = OVERLAP_N
               IF(FCHAIN_FN_MAG(FCHAINC)>FC_MAG_MAX) THEN
                  FC_MAX_INDEX = FCHAINC  ! Keep track of max magnitude for pair LL-I
                  FC_MAG_MAX = FCHAIN_FN_MAG(FCHAINC)
               ENDIF

            ENDIF


!Rolling friction
! JFD: Here the rolling friction coefficient is non-dimensional. This is the equivalent to
! the Zhou paper (Physica A 269 (1999) 536-553) 's definition (Model A) divided by the
! particle diameter. Therefore here mew_r is multiplied by the diameter.
            if(MEW_R>ZERO) then
               ! Particle LL
               OMEGA_MAG = sqrt(dot_product(OMEGA_NEW(LL,:),OMEGA_NEW(LL,:)))
               if(OMEGA_MAG > ZERO) then
                  RFTOW = -MEW_R*(2.0*DES_RADIUS(LL))*KN_DES * OVERLAP_N*OMEGA_NEW(LL,:)/OMEGA_MAG
                  TOW(LL,:) = TOW(LL,:) + RFTOW
               endif
               !particle i
               OMEGA_MAG = sqrt(dot_product(OMEGA_NEW(I,:),OMEGA_NEW(I,:)))
               if(OMEGA_MAG > ZERO) then
                  RFTOW = -MEW_R*(2.0*DES_RADIUS(I))*KN_DES * OVERLAP_N*OMEGA_NEW(I,:)/OMEGA_MAG
                  TOW(I,:) = TOW(I,:) + RFTOW
               endif
            endif

         ENDDO ! Neighbor CC i.e. particle I loop
         ! Update DES_COL_FORCE
         DES_COL_FORCE(LL,:) = FC(LL,:)

         IF(WRITE_FORCE_CHAIN.AND.FC_MAX_INDEX>0) FCHAIN_FCMAX(FC_MAX_INDEX) = .TRUE. ! Record where max magnitude occurs (for a given particle LL

      ENDDO ! Particle LL loop
!$omp end do

!$omp end parallel

! just for post-processing mag. of cohesive forces on each particle
      IF(USE_COHESION .AND. VAN_DER_WAALS .AND. GRAV_MAG > ZERO) THEN
         PostCohesive(:) = PostCohesive(:)/GRAV_MAG
      ENDIF

      RETURN

      contains

        include 'functions.inc'

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: print_excess_overlap                                    !
!                                                                      !
!  Purpose: Print overlap warning messages.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PRINT_EXCESS_OVERLAP

      use error_manager

      IF(OVERLAP_N > flag_overlap*DES_RADIUS(LL) .OR.                  &
         OVERLAP_N > flag_overlap*DES_RADIUS(I)) THEN

         WRITE(ERR_MSG,1000) trim(iVAL(LL)), trim(iVAL(I)), S_TIME,    &
            DES_RADIUS(LL), DES_RADIUS(I), OVERLAP_N

         CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1000 FORMAT('WARNING: Excessive overplay detected between ',          &
         'particles ',A,' and ',/A,' at time ',g11.4,'.',/             &
         'RADII:  ',g11.4,' and ',g11.4,4x,'OVERLAP: ',g11.4)

      END SUBROUTINE PRINT_EXCESS_OVERLAP

    END SUBROUTINE CALC_FORCE_DEM

 END MODULE CALC_FORCE_DEM_MOD
