#include "error.inc"

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_COLLISION_WALL                                    C
!  Author: Rahul Garg                               Date: 1-Dec-2013   C
!                                                                      C
!  Purpose: subroutines for particle-wall collisions when cutcell is   C
!           used. Also contains rehack of routines for cfslide and     C
!           cffctow which might be different from the stand alone      C
!           routines. Eventually all the DEM routines will be          C
!           consolidated.                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
MODULE CALC_COLLISION_WALL

   USE bc
   USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
   USE constant, only: pi
   USE cutcell, only: area_cut, blocked_cell_at, cartesian_grid
   USE des_thermo, only: flpc, calc_cond_des, des_t_s, q_source
   USE des_thermo_cond, only: des_conduction_wall, des_qw_cond
   USE discretelement
   USE error_manager
   USE exit, only: mfix_exit
   USE functions, only: funijk, fluid_at, is_normal
   USE geometry, only: do_k, imax1, imax2, imin1, imin2
   USE geometry, only: dx, dy, dz
   USE geometry, only: kmax1, kmax2, kmin1, kmin2, zlength
   USE geometry, only: no_k, jmax1, jmax2, jmin1, jmin2
   USE param, only: dimension_i, dimension_j, dimension_k
   USE param1, only: half, one, small_number, undefined
   USE particles_in_cell_mod, only: particles_in_cell
   USE physprop, only: k_g0, k_s0
   USE run, only: units, time
   USE stl
   USE stl_dbg_des, only: write_stls_this_dg
   USE stl_dbg_des, only: write_this_stl
   USE stl_functions_des, only: closestptpointtriangle

   use sq_contact_wall
   use SQ_EQUIVALENT_RADIUS_MOD
   use SQ_PROPERTIES_MOD

   PRIVATE
   PUBLIC :: CALC_DEM_FORCE_WITH_WALL_STL,&
      & CALC_DEM_THERMO_WITH_WALL_STL

CONTAINS

   SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL

      Implicit none

      INTEGER :: LL
      INTEGER :: NF
      DOUBLE PRECISION ::OVERLAP_N, SQRT_OVERLAP

      DOUBLE PRECISION :: V_REL_TRANS_NORM, DISTSQ, RADSQ, CLOSEST_PT(DIMN)
! local normal and tangential forces
      DOUBLE PRECISION :: NORMAL(DIMN), VREL_T(DIMN), DIST(DIMN), DISTMOD
      DOUBLE PRECISION, DIMENSION(DIMN) :: FTAN, FNORM, OVERLAP_T

      LOGICAL :: DES_LOC_DEBUG
      INTEGER :: CELL_ID, cell_count
      INTEGER :: PHASELL

      DOUBLE PRECISION :: TANGENT(DIMN)
      DOUBLE PRECISION :: FTMD, FNMD
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W

      double precision :: MAG_OVERLAP_T

      double precision :: line_t
! flag to tell if the orthogonal projection of sphere center to
! extended plane detects an overlap

      DOUBLE PRECISION :: MAX_DISTSQ, DISTAPART, FORCE_COH, R_LM
      INTEGER :: MAX_NF, axis
      DOUBLE PRECISION, DIMENSION(3) :: PARTICLE_MIN, PARTICLE_MAX, POS_TMP
      DOUBLE PRECISION, DIMENSION(3) :: NORM_PREVIOUS_FACET_COH
      DOUBLE PRECISION, DIMENSION(3) :: COH_FORCE_PREVIOUS_FACET
      DOUBLE PRECISION :: COHMAG_PREVIOUS_FACET
      DOUBLE PRECISION :: NORMTEST

! Flag to keep only cohesion force with one STL facet
      LOGICAL :: COHESION_WITH_STL_FACET
! Flag to distinguish point or edge intersection
      logical :: point_or_edge_int, moving_wall
! Superquadric DEM model
! Shape exponents of superquadric surface
      DOUBLE PRECISION :: m,n
! Euler parameters of superquadric surface, also called quaternions
      DOUBLE PRECISION :: euler0, euler1, euler2, euler3
! P is the location of superquadric center, axi is the sem-axes
      DOUBLE PRECISION :: p(3), axi(3),X(3)
      DOUBLE PRECISION :: contact_point_global(3),contact_point_wall_global(3)
      DOUBLE PRECISION :: norm_wall(3),vertex_wall(3),norm_wallsq,DIST2(3)
      DOUBLE PRECISION :: DIST3(3),dist_t
      integer :: contact_test_wall
! equivalent rad of local contact curvature
      DOUBLE PRECISION :: eq_r,r_eff_bs,r_eff_sq
      DOUBLE PRECISION :: HERT_KwN_SUPER(3),HERT_KwT_SUPER(3)
      DOUBLE PRECISION :: DES_ETAN_WALL_SUPER(3),DES_ETAT_WALL_SUPER(3)
      INTEGER :: shape_ll
      logical, parameter :: ldebug = .false.
! Rolling friction
      DOUBLE PRECISION :: RFTOW(3), omega_mag

      DES_LOC_DEBUG = .false. ;      DEBUG_DES = .false.
      FOCUS_PARTICLE = -1

! change the coefficient of restitution--yuxing
! Store the particle-pair imformation last time-step
      DOUBLE PRECISION :: DES_POS_L_OLD(3), DES_VEL_L_OLD(3), DES_OMEGA_L_OLD(3)
      DOUBLE PRECISION :: DIST_OLD(3), NORMAL_OLD(3), DISTSQ_OLD, V_REL_TRANS_NORM_OLD, VREL_T_OLD

!$omp parallel default(none) private(LL,MAG_OVERLAP_T,                    &
!$omp    cell_id,radsq,particle_max,particle_min,tangent,                 &
!$omp    axis,nf,closest_pt,dist,r_lm,distapart,force_coh,distsq,         &
!$omp    line_t,max_distsq,max_nf,normal,distmod,overlap_n,VREL_T,        &
!$omp    v_rel_trans_norm,phaseLL,sqrt_overlap,kn_des_w,kt_des_w,         &
!$omp    cohesion_with_stl_facet, point_or_edge_int, des_usr_var,         &
!$omp    norm_previous_facet_coh, cohmag_previous_facet,                  &
!$omp    coh_force_previous_facet, normtest,                              &
!$omp    etan_des_w,etat_des_w,fnorm,overlap_t,ftan,ftmd,fnmd,pos_tmp,    &
!$omp    mew_rw,rftow, omega_new,omega_mag,euler0,euler1,euler2,euler3,   &
!$omp    p,axi,m,n,norm_wall,norm_wallsq,vertex_wall)                     &
!$omp shared(max_pip,focus_particle,debug_des,                            &
!$omp    pijk,dg_pijk,des_pos_new,                                        &
!$omp    des_radius,facets_at_dg,vertex,                                  &
!$omp    ignore_edge_intersection, wall_surface_energy, cos_stl_nb_angle, &
!$omp    hert_kwn,hert_kwt,kn_w,kt_w,des_coll_model_enum,mew_w,tow,       &
!$omp    des_etan_wall,des_etat_wall,dtsolid,fc,norm_face,                &
!$omp    wall_collision_facet_id,wall_collision_PFT,use_cohesion,         &
!$omp    van_der_waals,wall_hamaker_constant,wall_vdw_outer_cutoff,       &
!$omp    wall_vdw_inner_cutoff,asperities,surface_energy,                 &
!$omp    super_q, super_mn, super_r,shape_ll, SuperDEM,                   &
!$omp    contact_test_wall,contact_point_wall_global,contact_point_global,&
!$omp    dist_t,eq_r,r_eff_sq,r_eff_bs,dist3,                             &
!$omp    des_etan_wall_super,des_etat_wall_super,                         &
!$omp    hert_kwn_super,hert_kwt_super)
!$omp do
      DO LL = 1, MAX_PIP

         IF(LL.EQ.FOCUS_PARTICLE) DEBUG_DES = .TRUE.

! skipping non-existent particles or ghost particles
! make sure the particle is not classified as a new 'entering' particle
! or is already marked as a potential exiting particle
         IF(.NOT.IS_NORMAL(LL)) CYCLE

         CELL_ID = DG_PIJK(LL)

! If no neighboring facet in the surrounding 26 cells, then exit
         IF(facets_at_dg(CELL_ID)%COUNT < 1) THEN
            WALL_COLLISION_FACET_ID(:,LL) = -1
            WALL_COLLISION_PFT(:,:,LL) = 0.0d0
            CYCLE
         ENDIF

         if (SuperDEM) then
! Quaternions for LL particle
             euler0= super_q(ll,1)
             euler1= super_q(ll,2)
             euler2= super_q(ll,3)
             euler3= super_q(ll,4)
! position of LL particle
             p(1)=DES_POS_NEW(LL,1)
             p(2)=DES_POS_NEW(LL,2)
             p(3)=DES_POS_NEW(LL,3)
! Semi-axis and roundness of LL particle
             axi(1)=super_r(ll,1)
             axi(2)=super_r(ll,2)
             axi(3)=super_r(ll,3)
             m=super_mn(ll,1)
             n=super_mn(ll,2)

! Check the shape of ll particle
             call Sq_shapes(axi,m,n,shape_ll)
         endif

! Check particle LL for wall contacts
         RADSQ = DES_RADIUS(LL)*DES_RADIUS(LL)

         particle_max(:) = des_pos_new( LL,:) + des_radius(LL)
         particle_min(:) = des_pos_new( LL,:) - des_radius(LL)

         COHESION_WITH_STL_FACET  = .FALSE.

         DO CELL_COUNT = 1, facets_at_dg(cell_id)%count

            axis = facets_at_dg(cell_id)%dir(cell_count)

            NF = facets_at_dg(cell_id)%id(cell_count)

! Compute particle-particle VDW cohesive short-range forces
            IF(USE_COHESION .AND. VAN_DER_WAALS) THEN

               CALL ClosestPtPointTriangle(DES_POS_NEW(LL,:),          &
                  VERTEX(:,:,NF), CLOSEST_PT(:),point_or_edge_int)

! To avoid counting cohesion forces multiple times when several
! triangles define a plane (or near planar) surface, edge or point
! interaction is ignored.
              IF(.NOT.(point_or_edge_int.AND.IGNORE_EDGE_INTERSECTION(NF))) THEN

                  DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(LL,:)
                  DISTSQ = DOT_PRODUCT(DIST, DIST)
                  R_LM = 1*DES_RADIUS(LL)

                  IF(DISTSQ < (R_LM+WALL_VDW_OUTER_CUTOFF)**2) THEN
                     IF(DISTSQ > (WALL_VDW_INNER_CUTOFF+R_LM)**2) THEN
                        DistApart = (SQRT(DISTSQ)-R_LM)
                        FORCE_COH = WALL_HAMAKER_CONSTANT*DES_RADIUS(LL) /&
                           (6.0d0*DistApart**2)*(Asperities/(Asperities +  &
                           DES_RADIUS(LL)) + ONE/(ONE+Asperities/         &
                           DistApart)**2)
                     ELSE
                        FORCE_COH = 4.0d0*PI*WALL_SURFACE_ENERGY*DES_RADIUS(LL)* &
                           (Asperities/(Asperities+DES_RADIUS(LL)) + ONE/ &
                           (ONE+Asperities/WALL_VDW_INNER_CUTOFF)**2 )
                     ENDIF


! The logic below is also aimed at avoiding counting the contribution of
! a wall multiple times. Once the contribution of a facet is accounted
! for, the contribution of the next facet is either:
! - ignored if the two normal vectors are close to each other (angle below STL_NB_ANGLE defined in
!   stl_mod.f), and its cohesion force magnitude is smaller (because only the closest
!   point should be kept)
! - kept, if the normals are close to each other and the cohesion force
!   magnitude is larger. In that case the force from the previous facet is
!   removed.
! If the facet is not part tof the same plane (say near a corner), the
! its contribution is added
                     IF(.NOT.COHESION_WITH_STL_FACET) THEN ! Record contribution of the first active facet
                       COHESION_WITH_STL_FACET = .TRUE.
                       NORM_PREVIOUS_FACET_COH = NORM_FACE(:,NF)
                       COHMAG_PREVIOUS_FACET = FORCE_COH/SQRT(DISTSQ)
                       COH_FORCE_PREVIOUS_FACET = DIST(:)*FORCE_COH/SQRT(DISTSQ)
                     ELSE
                        NORMTEST = DOT_PRODUCT(NORM_PREVIOUS_FACET_COH,NORM_FACE(:,NF))
                        IF(NORMTEST>COS_STL_NB_ANGLE) THEN
                           IF((FORCE_COH/SQRT(DISTSQ))>COHMAG_PREVIOUS_FACET) THEN
                              FC(LL,:) = FC(LL,:) - COH_FORCE_PREVIOUS_FACET ! Remove contribution from previous facet
                              COHESION_WITH_STL_FACET = .TRUE. ! and record contribution of current facet
                              NORM_PREVIOUS_FACET_COH = NORM_FACE(:,NF)
                              COHMAG_PREVIOUS_FACET = FORCE_COH/SQRT(DISTSQ)
                              COH_FORCE_PREVIOUS_FACET = DIST(:)*FORCE_COH/SQRT(DISTSQ)
                           ELSE
                              COHESION_WITH_STL_FACET = .FALSE.
                           ENDIF
                        ELSE  ! Keep adding contribution of current facet because it is not is he same plane (e.g, near a corner)
                           COHESION_WITH_STL_FACET = .TRUE.
                           NORM_PREVIOUS_FACET_COH = NORM_FACE(:,NF)
                           COHMAG_PREVIOUS_FACET = FORCE_COH/SQRT(DISTSQ)
                           COH_FORCE_PREVIOUS_FACET = DIST(:)*FORCE_COH/SQRT(DISTSQ)
                        ENDIF
                     ENDIF

                     IF(COHESION_WITH_STL_FACET) FC(LL,:) = FC(LL,:) + DIST(:)*FORCE_COH/SQRT(DISTSQ)


                  ENDIF
               ENDIF
            ENDIF

            if (facets_at_dg(cell_id)%min(cell_count) >    &
               particle_max(axis)) then
               call remove_collision(LL, nf, wall_collision_facet_id)
               cycle
            endif

            if (facets_at_dg(cell_id)%max(cell_count) <    &
               particle_min(axis)) then
               call remove_collision(LL, nf, wall_collision_facet_id)
               cycle
            endif

! Checking all the facets is time consuming due to the expensive
! separating axis test. Remove this facet from contention based on
! a simple orthogonal projection test.

! Parametrize a line as p = p_0 + t normal and intersect with the
! triangular plane. If t>0, then point is on the non-fluid side of
! the plane. If the plane normal is assumed to point toward the fluid.

! -undefined, because non zero values will imply the sphere center
! is on the non-fluid side of the plane. Since the testing
! is with extended plane, this could very well happen even
! when the particle is well inside the domain (assuming the plane
! normal points toward the fluid). See the pic below. So check
! only when line_t is negative

!                 \   Solid  /
!                  \  Side  /
!                   \      /
!                    \    /
! Wall 1, fluid side  \  /  Wall 2, fluid side
!                      \/
!                        o particle
!
! line_t will be positive for wall 1 (incorrectly indicating center
! is outside the domain) and line_t will be negative for wall 2.
!
! Therefore, only stick with this test when line_t is negative and let
! the separating axis test take care of the other cases.

! Since this is for checking static config, line's direction is the
! same as plane's normal. For moving particles, the line's normal will
! be along the point joining new and old positions.

            line_t = DOT_PRODUCT(VERTEX(1,:,NF) - des_pos_new(LL,:),&
               NORM_FACE(:,NF))

! k - rad >= tol_orth, where k = -line_t, then orthogonal
! projection is false. Substituting for k
! => line_t + rad <= -tol_orth
! choosing tol_orth = 0.01% of des_radius = 0.0001*des_radius

! Orthogonal projection will detect false positives even
! when the particle does not overlap the triangle.
! However, if the orthogonal projection shows no overlap, then
! that is a big fat negative and overlaps are not possible.

            if((line_t.le.-1.0001d0*des_radius(LL))) then  ! no overlap
                call remove_collision(LL,nf,wall_collision_facet_id)
                CYCLE
            ENDIF

! SuperDEM
            IF(SuperDEM .and. shape_ll .gt. 0) THEN
! Facet wall normal point to the fluid domain side
                norm_wall(:)=NORM_FACE(:,NF)
                norm_wallsq=(norm_wall(1))**2+(norm_wall(2))**2+(norm_wall(3))**2
                norm_wall(:)=norm_wall(:)/dsqrt(norm_wallsq)
                vertex_wall(:)=VERTEX(1,:,NF)

                call SQ_CONTACT_WITH_STL(p,axi,m,n,&
                    euler0,euler1, euler2, euler3,NORM_wall,&
                    vertex_wall,contact_point_global,&
                    contact_point_wall_global,contact_test_wall)

                pos_tmp(:) = contact_point_global(:)
                CALL ClosestPtPointTriangle(pos_tmp,             &
                    VERTEX(:,:,NF), CLOSEST_PT(:),point_or_edge_int)

                DIST(:) = contact_point_global(:)-contact_point_wall_global(:)
               dist_t=dot_product(dist(:), norm_wall(:))

               DISTSQ = dist_t**2.0


               if (contact_test_wall > 0) then  !not contact

                   call remove_collision(LL,nf,wall_collision_facet_id)
                   CYCLE
               endif

           ELSE ! SuperDEM

               pos_tmp = DES_POS_NEW(LL,:)
               CALL ClosestPtPointTriangle(pos_tmp,             &
                   VERTEX(:,:,NF), CLOSEST_PT(:),point_or_edge_int)

               DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(LL,:)
               DISTSQ = DOT_PRODUCT(DIST, DIST)

               IF(DISTSQ .GE. RADSQ - SMALL_NUMBER) THEN !No overlap exists
                   call remove_collision(LL,nf,wall_collision_facet_id)
                   CYCLE
               ENDIF

           ENDIF !SuperDEM
! If angle between two facets is small (below STL_NB_ANGLE defined in
! stl_mod.f), the collision with edge is ignored. This is done to avoid
! unwanted collision with all trianle edges making up a large plane or
! surface with low curvature
           IF(point_or_edge_int.AND.IGNORE_EDGE_INTERSECTION(NF)) THEN
               CYCLE
           ENDIF

           MAX_DISTSQ = DISTSQ
           MAX_NF = NF

! Assign the collision normal based on the facet with the
! largest overlap.

! Facet's normal is correct normal only when the intersection is with
! the face. When the intersection is with edge or vertex, then the
! normal is based on closest pt and sphere center. The definition above
! of the normal is generic enough to account for differences between
! vertex, edge, and facet.

! Calculate the particle/wall overlap.
           DISTMOD = SQRT(MAX_DISTSQ)
           if (SuperDEM .and. shape_ll .gt. 0) then
               OVERLAP_N = DISTMOD
               NORMAL(:)= - norm_wall(:)
           else
               NORMAL(:) = DIST(:)/sqrt(DISTSQ)
               OVERLAP_N = DES_RADIUS(LL) - DISTMOD
           endif

! Calculate the spring model parameters.
           phaseLL = PIJK(LL,5)

! Calculate the translational relative velocity
           CALL CFRELVEL_WALL_SuperDEM(LL,pos_tmp, p,V_REL_TRANS_NORM,&
               VREL_T,NORMAL, DISTMOD,shape_ll)
            
! Hertz vs linear spring-dashpot contact model
           IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               sqrt_overlap = SQRT(OVERLAP_N)

               IF(SuperDEM .and. shape_ll .gt. 0) THEN
! Equivalent radius of particle A at contact point
                   call EQUIVALENT_RADIUS(p,contact_point_global,axi,m,n,&
                       euler0,euler1,euler2,euler3, &
                       eq_r)
                   !
! Update the value of R_EFF for superquadric particle
                   R_EFF_sq = eq_r
                   R_EFF_bs=DES_RADIUS(LL)

! change the coefficient of restitution--yuxing
         ! calculate the particle-pair relative velocity last timestep: V_REL_TRANS_NORM_OLD, NORM_OLD
                   DES_POS_L_OLD = DES_POS_OLD_Eu(L,:)
                   DES_VEL_L_OLD = DES_VEL_OLD_Eu(L,:)
                   DES_OMEGA_L_OLD = OMEGA_OLD_Eu(L,:)

                   DIST_OLD(:) = CLOSEST_PT(:) - DES_POS_L_OLD
                   DISTSQ_OLD = DOT_PRODUCT(DIST_OLD, DIST_OLD)
                   NORMAL_OLD(:) = DIST_OLD(:)/sqrt(DISTSQ_OLD)

                   call Cal_Rel_Vel_Wall(DES_VEL_L_OLD, DES_OMEGA_L_OLD,  NORMAL_OLD, DISTSQ_OLD, V_REL_TRANS_NORM_OLD, VREL_T_OLD)
                   
         ! change the collosion paremeters based on velocity change: V_REL_TRANS_NORM_OLD, NORM_OLD, V_REL_TRANS_NORM, NORMAL(:)
                   if (( V_REL_TRANS_NORM_OLD .ge. 0.0) .and. (abs(V_REL_TRANS_NORM) .gt. abs(V_REL_TRANS_NORM_OLD)))  then
                     if ( V_REL_TRANS_NORM * V_REL_TRANS_NORM_OLD .ge. 0.0 ) then
                        IF(CALL_USR) then
                           ! CALL USR1_DES(abs(V_REL_TRANS_NORM),L,I)
                           call interpolate_keyframe_data(Current_V)
                           EN_Var(L) = kf_data(1)%var_current(1)
                           call change_LSD_parameters_Wall(L,I)
                        end if
                     end if
                  end if
      
                  if (( V_REL_TRANS_NORM_OLD .lt. 0.0) .and. (abs(V_REL_TRANS_NORM) .lt. abs(V_REL_TRANS_NORM_OLD)))  then
                     IF(CALL_USR) then
                        call interpolate_keyframe_data(Current_V)
                        EN_Var(L) = kf_data(1)%var_current(1)
                        call change_LSD_parameters_Wall(L,I)
                     end if
                  end if      

! Update some parameters for superquadric particle

                   DES_ETAN_WALL_SUPER(phaseLL)=DES_ETAN_WALL(phaseLL)/&
                       SQRT(HERT_KWN(phaseLL))
                   DES_ETAT_WALL_SUPER(phaseLL)= DES_ETAT_WALL(phaseLL)/&
                       SQRT(HERT_KWT(phaseLL))

                   hert_kWn_SUPER(phaseLL)=hert_kWn(phaseLL)/&
                       SQRT(R_eff_bs)*SQRT(R_eff_sq)
                   hert_kWt_SUPER(phaseLL)=hert_kWt(phaseLL)/&
                       SQRT(R_eff_bs)*SQRT(R_eff_sq)
                   DES_ETAN_wall_SUPER(phaseLL)=DES_ETAN_wall_SUPER(phaseLL)*&
                       SQRT(HERT_KwN_SUPER(phaseLL))
                   DES_ETAT_wall_SUPER(phaseLL)=DES_ETAT_wall_SUPER(phaseLL)*&
                       SQRT(HERT_KwT_SUPER(phaseLL))

                   KN_DES_W = hert_kWn_super(phaseLL)*sqrt_overlap
                   KT_DES_W = hert_kWt_super(phaseLL)*sqrt_overlap
                   sqrt_overlap = SQRT(sqrt_overlap)
                   ETAN_DES_W = DES_ETAN_WALL_super(phaseLL)*sqrt_overlap
                   ETAT_DES_W = DES_ETAT_WALL_super(phaseLL)*sqrt_overlap

               ELSE  ! for sphere, hertz

                   KN_DES_W = hert_kwn(phaseLL)*sqrt_overlap
                   KT_DES_W = hert_kwt(phaseLL)*sqrt_overlap
                   sqrt_overlap = SQRT(sqrt_overlap)
                   ETAN_DES_W = DES_ETAN_WALL(phaseLL)*sqrt_overlap
                   ETAT_DES_W = DES_ETAT_WALL(phaseLL)*sqrt_overlap
               ENDIF

           ELSE ! for sphere, lsd
               KN_DES_W = KN_W
               KT_DES_W = KT_W
               ETAN_DES_W = DES_ETAN_WALL(phaseLL)
               ETAT_DES_W = DES_ETAT_WALL(phaseLL)
           ENDIF

           FNORM(:) = -(KN_DES_W * OVERLAP_N * NORMAL(:) + &
               ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:))

           if (ldebug) then
               WRITE(*,*) 'norm force and v_rel_trans_norm'
               WRITE(*,'(3(f18.5))')   FNORM(:)
               WRITE(*,'(3(f18.5))') V_REL_TRANS_NORM,overlap_n
           endif
! Calculate the tangential displacement.
           OVERLAP_T(:) = DTSOLID*VREL_T(:) + GET_COLLISION(LL,       &
               NF, WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)
           MAG_OVERLAP_T = sqrt(DOT_PRODUCT(OVERLAP_T, OVERLAP_T))
           if (ldebug) then
               WRITE(*,*) 'overlaps_t and mag_overlap_t'
               WRITE(*,'(3(f18.5))')   overlap_t
               WRITE(*,'(3(f18.5))') mag_overlap_t
           endif

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall.
           IF(MAG_OVERLAP_T > 0.0) THEN
! Calculate the tangential contact force.
               FTAN = -KT_DES_W*OVERLAP_T - ETAT_DES_W*VREL_T
               FTMD = sqrt(dot_product(FTAN,FTAN))
! Max force before the on set of frictional slip.
               FNMD = MEW_W*sqrt(dot_product(FNORM,FNORM))
! Frictional slip
               IF(FTMD > FNMD) THEN
! Direction of tangential force.
                   TANGENT = OVERLAP_T/MAG_OVERLAP_T
                   FTAN = -FNMD * TANGENT
                   OVERLAP_T = (FNMD/KT_DES_W) * TANGENT
               ENDIF
           ELSE
               FTAN = 0.0
           ENDIF

! Save the tangential displacement.
           CALL UPDATE_COLLISION(OVERLAP_T, LL, NF,                   &
               WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)

! Add the collision force to the total forces acting on the particle.
           FC(LL,:) = FC(LL,:) + FNORM(:) + FTAN(:) 

! Add the torque force to the total torque acting on the particle.
           IF (SuperDEM .and. shape_ll .gt. 0) THEN
               DIST3(:) = contact_point_wall_global(:) - DES_POS_NEW(LL,:)
               TOW(LL,:) = TOW(LL,:) + DES_CROSSPRDCT(DIST3,(FTAN+FNORM))
           ELSE

               TOW(LL,:) = TOW(LL,:) + DISTMOD*DES_CROSSPRDCT(NORMAL,FTAN)
           ENDIF

!Rolling friction
! JFD: Here the rolling friction coefficient is non-dimensional. This is the equivalent to
! the definition from Zhou (Physica A 269 (1999) 536-553) divided by the
! particle diameter. Therefore here mew_rw is multiplied by the diameter.
            if(MEW_RW>ZERO) then
               OMEGA_MAG = sqrt(dot_product(OMEGA_NEW(LL,:),OMEGA_NEW(LL,:)))
               if(OMEGA_MAG > ZERO) then
                 RFTOW = -MEW_RW*(2.0*DES_RADIUS(LL))*KN_DES_W * OVERLAP_N*OMEGA_NEW(LL,:)/OMEGA_MAG
                 TOW(LL,:) = TOW(LL,:) + RFTOW
               endif
            endif

       ENDDO

   ENDDO
!$omp end do
!$omp end parallel

   RETURN

contains

!......................................................................!
!  Function: GET_COLLISION                                             !
!                                                                      !
!  Purpose: Return the integrated (t0->t) tangential displacement.     !
!......................................................................!
     FUNCTION GET_COLLISION(LLL,FACET_ID,WALL_COLLISION_FACET_ID,     &
         WALL_COLLISION_PFT)

          IMPLICIT NONE

      DOUBLE PRECISION :: GET_COLLISION(DIMN)
      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, allocatable, INTENT(INOUT) :: WALL_COLLISION_FACET_ID(:,:)
      DOUBLE PRECISION, allocatable, INTENT(INOUT) :: WALL_COLLISION_PFT(:,:,:)
      INTEGER :: CC, FREE_INDEX, LC, dgIJK


      free_index = -1

      do cc = 1, COLLISION_ARRAY_MAX
         if (facet_id == wall_collision_facet_id(cc,LLL)) then
            get_collision(:) = wall_collision_PFT(:,cc,LLL)
            return
         else if (-1 == wall_collision_facet_id(cc,LLL)) then
            free_index = cc
         endif
      enddo

! Overwrite old data. This is needed because a particle moving from
! one dg cell to another may no longer 'see' an STL before it moved
! out of contact range. Therefore, the 'remove_collision' function
! does not get called to cleanup the stale data.
      if(-1 == free_index) then
         dgIJK=DG_PIJK(LLL)
         cc_lp: do cc=1, COLLISION_ARRAY_MAX
            do lc=1, facets_at_dg(dgIJK)%count
               if(wall_collision_facet_id(cc,LLL) == &
                  facets_at_dg(dgIJK)%id(LC))  cycle cc_lp
            enddo
            free_index = cc
            exit cc_lp
         enddo cc_lp
      endif

! Last resort... grow the collision array
      if(-1 == free_index) then
         free_index=COLLISION_ARRAY_MAX+1
         COLLISION_ARRAY_MAX = 2*COLLISION_ARRAY_MAX
         CALL GROW_WALL_COLLISION(COLLISION_ARRAY_MAX)
      endif

      wall_collision_facet_id(free_index,LLL) = facet_id
      wall_collision_PFT(:,free_index,LLL) = ZERO
      get_collision(:) = wall_collision_PFT(:,free_index,LLL)
      return

      END FUNCTION GET_COLLISION


!......................................................................!
!  Subroutine: GROW_WALL_COLLISION                                     !
!                                                                      !
!  Purpose: Return the integrated (t0->t) tangential displacement.     !
!......................................................................!
      SUBROUTINE GROW_WALL_COLLISION(NEW_SIZE)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NEW_SIZE
      INTEGER :: lSIZE1, lSIZE2, lSIZE3
      INTEGER, ALLOCATABLE :: tmpI2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: tmpR3(:,:,:)

      lSIZE1 = size(wall_collision_facet_id,1)
      lSIZE2 = size(wall_collision_facet_id,2)

      allocate(tmpI2(NEW_SIZE, lSIZE2))
      tmpI2(1:lSIZE1,:) = WALL_COLLISION_FACET_ID(1:lSIZE1,:)
      call move_alloc(tmpI2, WALL_COLLISION_FACET_ID)
      WALL_COLLISION_FACET_ID(lSIZE1+1:NEW_SIZE,:) = -1

      lSIZE1 = size(wall_collision_pft,1)
      lSIZE2 = size(wall_collision_pft,2)
      lSIZE3 = size(wall_collision_pft,3)

      allocate(tmpR3(lSIZE1, NEW_SIZE, lSIZE3))
      tmpR3(:,1:lSIZE2,:) = WALL_COLLISION_PFT(:,1:lSIZE2,:)
      call move_alloc(tmpR3, WALL_COLLISION_PFT)

      RETURN
      END SUBROUTINE GROW_WALL_COLLISION




!......................................................................!
!  Function: UPDATE_COLLISION                                          !
!                                                                      !
!  Purpose: Update the integrated (t0->t) tangential displacement.     !
!......................................................................!
      SUBROUTINE UPDATE_COLLISION(PFT, LLL, FACET_ID,                  &
         WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)

      implicit none

      DOUBLE PRECISION, INTENT(IN) :: PFT(DIMN)
      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(IN) :: WALL_COLLISION_FACET_ID(:,:)
      DOUBLE PRECISION, INTENT(INOUT) :: WALL_COLLISION_PFT(:,:,:)
      INTEGER :: CC

      do cc = 1, COLLISION_ARRAY_MAX
         if (facet_id == wall_collision_facet_id(cc,LLL)) then
            wall_collision_PFT(:,cc,LLL) = PFT(:)
            return
         endif
      enddo

      WRITE(ERR_MSG, 1100)
      CALL LOG_ERROR()

 1100 FORMAT('Error: COLLISION_ARRAY_MAX too small. ')

      END SUBROUTINE UPDATE_COLLISION

!......................................................................!
!  Function: REMOVE_COLLISION                                          !
!                                                                      !
!  Purpose: Clear the integrated (t0->t) tangential displacement once  !
!  the collision is over (contact ended).                              !
!......................................................................!
      SUBROUTINE REMOVE_COLLISION(LLL,FACET_ID,WALL_COLLISION_FACET_ID)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(INOUT) :: WALL_COLLISION_FACET_ID(:,:)
      INTEGER :: CC

      DO CC = 1, COLLISION_ARRAY_MAX
         IF (FACET_ID == WALL_COLLISION_FACET_ID(CC,LLL)) THEN
            WALL_COLLISION_FACET_ID(CC,LLL) = -1
            RETURN
         ENDIF
      ENDDO

      END SUBROUTINE REMOVE_COLLISION

      END SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CFRELVEL_WALL                                           !
!                                                                      !
!  Purpose: Calculate the normal and tangential components of the      !
!  relative velocity between a particle and wall contact.              !
!                                                                      !
!  Comments: Only the magnitude of the normal component is returned    !
!  whereas the full tangential vector is returned.                     !
!----------------------------------------------------------------------!
      SUBROUTINE CFRELVEL_WALL(LL, VRN, VRT, NORM, DIST)

      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Particle index.
      INTEGER, INTENT(IN) :: LL
! Magnitude of the total relative translational velocity.
      DOUBLE PRECISION, INTENT(OUT):: VRN
! Total relative translational velocity (vector).
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(OUT):: VRT
! Unit normal from particle center to closest point on stl (wall)
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: NORM
! Distance between particle center and stl (wall).
      DOUBLE PRECISION, INTENT(IN) :: DIST

! Local variables
!---------------------------------------------------------------------//
! Additional relative translational motion due to rotation
      DOUBLE PRECISION, DIMENSION(DIMN) :: V_ROT
! Total relative velocity at contact point
      DOUBLE PRECISION, DIMENSION(DIMN) :: VRELTRANS

! Total relative velocity + rotational contribution
      V_ROT = DIST*OMEGA_NEW(LL,:)
      VRELTRANS(:) =  DES_VEL_NEW(LL,:) + DES_CROSSPRDCT(V_ROT, NORM)

! magnitude of normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! total relative translational slip velocity at the contact point
! Equation (8) in Tsuji et al. 1992
      VRT(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
      END SUBROUTINE CFRELVEL_WALL

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CFRELVEL_WALL_SuperDEM                                           !
!                                                                      !
!  Purpose: Calculate the normal and tangential components of the      !
!  relative velocity between a particle and wall contact.              !
!                                                                      !
!  Comments: Only the magnitude of the normal component is returned    !
!  whereas the full tangential vector is returned.                     !
!----------------------------------------------------------------------!
      SUBROUTINE CFRELVEL_WALL_SuperDEM(LL, X,POS_LL,VRN, VRT, NORM, DIST, shape_ll)
      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Particle index.
      INTEGER, INTENT(IN) :: LL
! Magnitude of the total relative translational velocity.
      DOUBLE PRECISION, INTENT(OUT):: VRN
! Total relative translational velocity (vector).
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(OUT):: VRT
! Unit normal from particle center to closest point on stl (wall)
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: NORM
! Distance between particle center and stl (wall).
      DOUBLE PRECISION, INTENT(IN) :: DIST
! momentum arms,Contact point
      DOUBLE PRECISION :: MOM_LL(3),X(3),POS_LL(3)
! Local variables
!---------------------------------------------------------------------//
! Additional relative translational motion due to rotation
      DOUBLE PRECISION, DIMENSION(DIMN) :: V_ROT
! Total relative velocity at contact point
      DOUBLE PRECISION, DIMENSION(DIMN) :: VRELTRANS
      INTEGER :: shape_ll

      if (SuperDEM .and. shape_ll .gt. 0) then
! Momentum of L particle  (center to contact point)
        MOM_LL(1)=X(1)-POS_LL(1)
        MOM_LL(2)=X(2)-POS_LL(2)
        MOM_LL(3)=X(3)-POS_LL(3)
        V_ROT(:) = DES_CROSSPRDCT(OMEGA_NEW(LL,:), MOM_LL(:))
      else
! Total relative velocity + rotational contribution
        V_ROT = DIST*OMEGA_NEW(LL,:)
        V_ROT = DES_CROSSPRDCT(V_ROT, NORM)
      endif

      VRELTRANS(:) =  DES_VEL_NEW(LL,:) + V_ROT(:)

! magnitude of normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! total relative translational slip velocity at the contact point
! Equation (8) in Tsuji et al. 1992
      VRT(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
      END SUBROUTINE CFRELVEL_WALL_SuperDEM


!----------------------------------------------------------------//
! SUBROUTINE: CALC_DEM_THERMO_WITH_WALL_STL
! By: Aaron M.
! Purpose: Compute heat transfer to particles due to boundaries
!----------------------------------------------------------------//

      SUBROUTINE CALC_DEM_THERMO_WITH_WALL_STL

      IMPLICIT NONE

      INTEGER :: LL             ! Loop index for particle
      INTEGER :: CELL_ID        ! Desgrid cell index
      INTEGER :: CELL_COUNT     ! Loop index for facets
      INTEGER :: IJK_FLUID      ! IJK index for fluid cell
      INTEGER :: I_FLUID, J_FLUID, K_FLUID
      INTEGER :: I_Facet, J_Facet, K_Facet, IJK_Facet
      INTEGER :: I1,J1,K1
      INTEGER :: phase_LL       ! Phase index for particle

      INTEGER, PARAMETER :: MAX_CONTACTS = 12
      INTEGER :: iFacet         ! loop index for facets
      INTEGER :: count_Facets   ! counter for number of facets particle contacts
      DOUBLE PRECISION, DIMENSION(3,MAX_CONTACTS) :: NORM_FAC_CONTACT
      LOGICAL :: USE_FACET
      LOGICAL :: DOMAIN_BDRY

      INTEGER :: NF             ! Facet ID
      INTEGER :: AXIS           ! Facet direction
      DOUBLE PRECISION :: RLENS_SQ ! lens radius squared
      DOUBLE PRECISION :: RLENS ! lens radius
      DOUBLE PRECISION :: RPART ! Particle radius

      ! vectors for minimum / maximum possible particle contacts
      DOUBLE PRECISION, DIMENSION(3) :: PARTPOS_MIN, PARTPOS_MAX
      DOUBLE PRECISION, DIMENSION(3) :: POS_TMP ! Temp vector
      DOUBLE PRECISION, DIMENSION(3) :: DIST, NORMAL

      DOUBLE PRECISION, DIMENSION(3) :: CLOSEST_PT

      DOUBLE PRECISION :: line_t  ! Normal projection for part/wall
      DOUBLE PRECISION :: DISTSQ ! Separation from particle and facet (squared)
      DOUBLE PRECISION :: PROJ

      INTEGER :: BC_ID ! BC ID
      INTEGER :: IBC   ! BC Loop Index

      DOUBLE PRECISION :: TWALL
      DOUBLE PRECISION :: K_gas
      DOUBLE PRECISION :: TPART
      DOUBLE PRECISION :: OVERLAP
      DOUBLE PRECISION :: QSWALL, AREA

      LOGICAL, SAVE :: OUTPUT_WARNING = .TRUE.

! Flag to distinguish point or edge intersection
      logical :: point_or_edge_int


    !  IF(.NOT.DES_CONTINUUM_COUPLED.or.DES_EXPLICITLY_COUPLED)THEN
    !     CALL PARTICLES_IN_CELL
    !  ENDIF
      DO LL = 1, MAX_PIP
         ! Skip non-existent particles or ghost particles
         IF (.NOT.IS_NORMAL(LL)) CYCLE
         PHASE_LL = PIJK(LL,5)
         IF(.NOT.CALC_COND_DES(PHASE_LL))CYCLE
         ! Get desgrid cell index
         CELL_ID = DG_PIJK(LL)

         ! Skip cells that do not have neighboring facet
         IF (facets_at_dg(CELL_ID)%COUNT <1) CYCLE

         ! Store lens radius of particle
         RLENS = (ONE+FLPC)*DES_RADIUS(LL)
         RLENS_SQ = RLENS*RLENS

         RPART = DES_RADIUS(LL)

         ! Compute max/min particle locations
         PARTPOS_MAX(:) = des_pos_new(LL,:) + RLENS
         PARTPOS_MIN(:) = des_pos_new(LL,:) - RLENS

         ! Get fluid cell
         I_FLUID = PIJK(LL,1)
         J_FLUID = PIJK(LL,2)
         K_FLUID = PIJK(LL,3)
         IJK_FLUID = PIJK(LL,4)
         TPART = DES_T_S(LL)

         ! Sometimes PIJK is not updated every solids timestep and
         ! therefore doesn't give the correct fluid cell that
         ! contains the particles.  If so, determine actual fluid
         ! cell that contains the particle
         IF(.NOT.DES_CONTINUUM_COUPLED.or.DES_EXPLICITLY_COUPLED)THEN
            IF(I_FLUID <= ISTART3 .OR. I_FLUID >= IEND3) THEN
               CALL PIC_SEARCH(I_FLUID, DES_POS_NEW(LL,1), XE,       &
               DIMENSION_I, IMIN2, IMAX2)
            ELSE
               IF((DES_POS_NEW(LL,1) >= XE(I_FLUID-1)) .AND.         &
               (DES_POS_NEW(LL,1) <  XE(I_FLUID))) THEN
               I_FLUID = I_FLUID
               ELSEIF((DES_POS_NEW(LL,1) >= XE(I_FLUID)) .AND.       &
                  (DES_POS_NEW(LL,1) < XE(I_FLUID+1))) THEN
                  I_FLUID = I_FLUID+1
               ELSEIF((DES_POS_NEW(LL,1) >= XE(I_FLUID-2)) .AND.     &
                  (DES_POS_NEW(LL,1) < XE(I_FLUID-1))) THEN
                  I_FLUID = I_FLUID-1
               ELSE
                  CALL PIC_SEARCH(I_FLUID, DES_POS_NEW(LL,1), XE,    &
                  DIMENSION_I, IMIN2, IMAX2)
               ENDIF
            ENDIF !(I)
            IF(J_FLUID <= JSTART3 .OR. J_FLUID >= JEND3) THEN
               CALL PIC_SEARCH(J_FLUID, DES_POS_NEW(LL,2), YN,       &
               DIMENSION_J, JMIN2, JMAX2)
            ELSE
               IF((DES_POS_NEW(LL,2) >= YN(J_FLUID-1)) .AND.         &
                  (DES_POS_NEW(LL,2) <  YN(J_FLUID))) THEN
                  J_FLUID = J_FLUID
               ELSEIF((DES_POS_NEW(LL,2) >= YN(J_FLUID)) .AND.       &
                  (DES_POS_NEW(LL,2) < YN(J_FLUID+1))) THEN
                  J_FLUID = J_FLUID+1
               ELSEIF((DES_POS_NEW(LL,2) >= YN(J_FLUID-2)) .AND.     &
                  (DES_POS_NEW(LL,2) < YN(J_FLUID-1))) THEN
                  J_FLUID = J_FLUID-1
               ELSE
                  CALL PIC_SEARCH(J_FLUID, DES_POS_NEW(LL,2), YN,    &
                  DIMENSION_J, JMIN2, JMAX2)
               ENDIF
            ENDIF !(J)

            IF(NO_K) THEN
               K_FLUID = 1
            ELSE
               IF(K_FLUID <= KSTART3 .OR. K_FLUID >= KEND3) THEN
               CALL PIC_SEARCH(K_FLUID, DES_POS_NEW(LL,3), ZT,         &
                  DIMENSION_K, KMIN2, KMAX2)
               ELSE
                  IF((DES_POS_NEW(LL,3) >= ZT(K_FLUID-1)) .AND.        &
                     (DES_POS_NEW(LL,3) < ZT(K_FLUID))) THEN
                     K_FLUID = K_FLUID
                  ELSEIF((DES_POS_NEW(LL,3) >= ZT(K_FLUID)) .AND.      &
                     (DES_POS_NEW(LL,3) < ZT(K_FLUID+1))) THEN
                     K_FLUID = K_FLUID+1
                  ELSEIF((DES_POS_NEW(LL,3) >= ZT(K_FLUID-2)) .AND.    &
                     (DES_POS_NEW(LL,3) < ZT(K_FLUID-1))) THEN
                     K_FLUID = K_FLUID-1
                  ELSE
                     CALL PIC_SEARCH(K_FLUID, DES_POS_NEW(LL,3), ZT,   &
                     DIMENSION_K, KMIN2, KMAX2)
                  ENDIF
               ENDIF
            ENDIF !(K)
            IJK_FLUID = PIJK(LL,4)
         ENDIF ! (NOT CONTINUUM_COUPLED OR EXPLICIT)


! Initialize counter for number of facets
         COUNT_FACETS = 0

         ! Loop over potential facets
         DO CELL_COUNT = 1, facets_at_dg(CELL_ID)%count
            ! Get direction (axis) and facet id (NF)
            axis = facets_at_dg(cell_id)%dir(cell_count)
            NF = facets_at_dg(cell_id)%id(cell_count)

            ! Check to see if facet is out of reach of particle
            if (facets_at_dg(cell_id)%min(cell_count) >    &
               partpos_max(axis)) then
               cycle
            endif
            if (facets_at_dg(cell_id)%max(cell_count) <    &
               partpos_min(axis)) then
               cycle
            endif

            ! Compute projection (normal distance) from wall to particle
            line_t = DOT_PRODUCT(VERTEX(1,:,NF) - des_pos_new(LL,:),&
            &        NORM_FACE(:,NF))

            ! If normal exceeds particle radius, particle is not in contact
            if((line_t.lt.-RLENS))CYCLE

            ! Compute closest point on facet
            POS_TMP(:) = DES_POS_NEW(LL,:)
            CALL ClosestPtPointTriangle(POS_TMP, VERTEX(:,:,NF), &
            &    CLOSEST_PT(:),point_or_edge_int)
            ! Compute position vector from particle to closest point on facet
            DIST(:) = CLOSEST_PT(:)-POS_TMP(:)
            DISTSQ = DOT_PRODUCT(DIST,DIST)

            ! Skip particles that are more than lens radius from facet
            IF(DISTSQ .GE. (RLENS_SQ-SMALL_NUMBER))CYCLE

            ! Only do heat transfer for particles that are directly above facet
            ! Heat transfer routines not generalized for edge contacts, but
            ! those should normally yield negligible heat transfer contributions

            ! Normalize distance vector and compare to facet normal
            NORMAL(:)=-DIST(:)/sqrt(DISTSQ)
            PROJ = sqrt(abs(DOT_PRODUCT(NORMAL, NORM_FACE(:,NF))))
            IF(ABS(ONE-PROJ).gt.1.0D-6)CYCLE

            ! Get overlap
            OVERLAP = RPART - SQRT(DISTSQ)

            ! Initialize area for facet (for post-proc. flux)
            AREA = ZERO

            ! Initialize BDRY FLAG
            DOMAIN_BDRY = .FALSE.
            ! Get BC_ID
            BC_ID = BC_ID_STL_FACE(NF)

            ! BC_ID is set to 0 in pre-proc if stl is a domain boundary
            IF(BC_ID.eq.0)then
               I1=I_FLUID
               J1=J_FLUID
               K1=K_FLUID
               DOMAIN_BDRY = .TRUE.

               ! Domain boundary, figure out real boundary ID
               IF(NORM_FACE(1,NF).ge.0.9999)THEN
                  ! WEST face
                  I1=IMIN2
                  IF(DO_K)THEN
                     AREA = DY(J_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DY(J_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(1,NF).le.-0.9999)THEN
                  ! EAST FACE
                  I1=IMAX2
                  IF(DO_K)THEN
                     AREA = DY(J_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DY(J_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(2,NF).ge.0.9999)THEN
                  ! SOUTH FACE
                  J1=JMIN2
                  IF(DO_K)THEN
                     AREA = DX(I_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DX(I_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(2,NF).le.-0.9999)THEN
                  ! NORTH FACE
                  J1=JMAX2
                  IF(DO_K)THEN
                     AREA = DX(I_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DX(I_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(3,NF).ge.0.9999)THEN
                  ! BOTTOM FACE
                  K1=KMIN2
                  AREA = DX(I_FLUID)*DY(J_FLUID)

               ELSEIF(NORM_FACE(3,NF).le.-0.9999)THEN
                  ! TOP FACE
                  K1=KMAX2
                  AREA = DX(I_FLUID)*DY(J_FLUID)
               ELSE
                  WRITE( *,*)'PROBLEM, COULD NOT FIND DOMAIN BOUNDARY'
                  WRITE(*,*)' In calc_thermo_des_wall_stl'
                  call mfix_exit(1)

               ENDIF


               ! Loop through defined BCs to see which one particle neighbors
               DO IBC = 1, DIMENSION_BC
                  IF(.NOT.BC_DEFINED(IBC))CYCLE
                  IF (I1.ge.BC_I_W(IBC).and.I1.le.BC_I_E(IBC).and.&
                      J1.ge.BC_J_S(IBC).and.J1.le.BC_J_N(IBC).and.&
                      K1.ge.BC_K_B(IBC).and.K1.le.BC_K_T(IBC))THEN
                      BC_ID = IBC
                      exit
                  ENDIF
               ENDDO
               IF(BC_ID.eq.0)then
                  IF(OUTPUT_WARNING)THEN
1111                  FORMAT("Warning: Could not find BC."/"Check input file to make sure domain boundaries are defined.",/ &
                          "DES_POS_NEW: ", 3(F12.5, 3X),/ &
                          "I,J,K: ", I4,I4,I4,          / &
                          "CLOSEST_PT: ", 3(F12.5, 3X), / &
                          "NORM_FACE: ",  3(F12.5, 3X), / &
                          'Suppressing further warnings.')
                      WRITE(ERR_MSG, 1111) (DES_POS_NEW(LL,IBC),IBC=1,3), &
                          I1, J1, K1,                                     &
                          (CLOSEST_PT(IBC),IBC=1,3),                      &
                          (NORM_FACE(IBC,NF),IBC=1,3)
                      OUTPUT_WARNING = .FALSE.
                      CALL LOG_WARNING()
                  ENDIF
                  CYCLE  !
               ENDIF

            ENDIF !Domain Boundary (facet ID was 0)

            IF (BC_TYPE_ENUM(BC_ID) == NO_SLIP_WALL .OR. &
               BC_TYPE_ENUM(BC_ID) == FREE_SLIP_WALL .OR. &
               BC_TYPE_ENUM(BC_ID) == PAR_SLIP_WALL .OR. &
               BC_TYPE_ENUM(BC_ID) == CG_NSW .OR. &
               BC_TYPE_ENUM(BC_ID) == CG_FSW .OR. &
               BC_TYPE_ENUM(BC_ID) == CG_PSW) THEN

               ! CHECK TO MAKE SURE FACET IS UNIQUE
               USE_FACET=.TRUE.
               DO IFACET=1,count_facets
                  ! DO CHECK BY ENSURING NORMAL VECTOR IS NEARLY PARALLEL
                  PROJ = sqrt(abs(DOT_PRODUCT(NORMAL, NORM_FAC_CONTACT(:,IFACET))))
                  IF(ABS(ONE-PROJ).lt.1.0D-6)THEN
                     USE_FACET=.FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF(.NOT.USE_FACET)CYCLE

               ! FACET IS UNIQUE
               count_facets=count_facets+1
               NORM_FAC_CONTACT(:,count_facets)=NORMAL(:)

! Do heat transfer
               ! GET WALL TEMPERATURE
               TWALL = BC_TW_S(BC_ID,phase_LL)

               ! GET GAS THERMAL CONDUCTIVITY
               if(k_g0.eq.UNDEFINED)then
                  ! Compute gas conductivity as is done in calc_k_g
                  ! But use average of particle and wall temperature to be gas temperature

                  K_Gas = 6.02D-5*SQRT(HALF*(TWALL+TPART)/300.D0) ! cal/(s.cm.K)
                  ! 1 cal = 4.183925D0 J
                  IF (UNITS == 'SI') K_Gas = 418.3925D0*K_Gas !J/s.m.K
               else
                  K_Gas=k_g0
               endif
               IF(TWALL.eq.UNDEFINED)CYCLE
               QSWALL = DES_CONDUCTION_WALL(LL, OVERLAP,K_s0(phase_LL), &
               &        K_s0(phase_LL),K_Gas,TWALL, TPART, RPART, &
               &        RLENS, phase_LL)


               Q_Source(LL) = Q_Source(LL)+QSWALL

               ! BELOW CODE IS ONLY NECESSARY FOR OUTPUTTING
               ! DATA.  Need to know fluid cell that contact
               ! point resides in so that wall flux can be
               ! output correctly.

               I_FACET = I_FLUID
               J_FACET = J_FLUID
               K_FACET = K_FLUID

               ! This checks to see if the contact point was NOT on a domain boundary

               IF(.NOT.DOMAIN_BDRY)THEN
                  IF(CLOSEST_PT(1) >= XE(I_FLUID))THEN
                     I_FACET = MIN(IMAX1, I_FLUID+1)
                  ELSEIF(CLOSEST_PT(1) < XE(I_FLUID-1))THEN
                     I_FACET = MAX(IMIN1, I_FLUID-1)
                  ENDIF

                  IF(CLOSEST_PT(2) >= YN(J_FLUID))THEN
                     J_FACET = MIN(JMAX1, J_FLUID+1)
                  ELSEIF(CLOSEST_PT(2) < YN(J_FLUID-1))THEN
                     J_FACET = MAX(JMIN1, J_FLUID-1)
                  ENDIF
                  IF(DO_K)THEN
                     IF(CLOSEST_PT(3) >= ZT(K_FLUID))THEN
                        K_FACET = MIN(KMAX1, K_FLUID+1)
                     ELSEIF(CLOSEST_PT(3) < ZT(K_FLUID-1))THEN
                        K_FACET = MAX(KMIN1, K_FLUID-1)
                     ENDIF
                  ENDIF
                  IJK_FACET=funijk(I_facet,J_facet,K_facet)
                  if (cartesian_grid) AREA=AREA_CUT(IJK_FACET)

               ! AREA is left undefined if contact was with cut-cell surface
                  IF (AREA.eq.ZERO)then
                     I_FACET = I_FLUID
                     J_FACET = J_FLUID
                     K_FACET = K_FLUID
                     IJK_FACET=funijk(I_facet,J_facet,K_facet)
                     IF(NORM_FACE(1,NF).ge.0.9999)THEN
                     ! WEST face
                        IF(DO_K)THEN
                           AREA = DY(J_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DY(J_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(1,NF).le.-0.9999)THEN
                     ! EAST FACE
                        IF(DO_K)THEN
                           AREA = DY(J_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DY(J_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(2,NF).ge.0.9999)THEN
                     ! SOUTH FACE
                        IF(DO_K)THEN
                           AREA = DX(I_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DX(I_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(2,NF).le.-0.9999)THEN
                     ! NORTH FACE
                        IF(DO_K)THEN
                           AREA = DX(I_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DX(I_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(3,NF).ge.0.9999)THEN
                     ! BOTTOM FACE
                        AREA = DX(I_FACET)*DY(J_FACET)
                     ELSEIF(NORM_FACE(3,NF).le.-0.9999)THEN
                     ! TOP FACE
                        AREA = DX(I_FACET)*DY(J_FACET)
                     ENDIF
                  ENDIF ! Area==0 (because cut-cell facet exists in non cut-cell)
               ENDIF ! NOT DOMAIN_BDRY
               IJK_FACET = FUNIJK(I_FACET,J_FACET,K_FACET)
               ! AM error check
               ! todo use proper error handling
               IF(cartesian_grid) THEN
                  IF(.NOT.FLUID_AT(IJK_FACET).AND. &
                     .NOT.BLOCKED_CELL_AT(IJK_FACET)) THEN
                     write(*,*)'ERROR: Cell containing facet is not a fluid', &
                     'or blocked cell'
                     write(*,*)FLUID_AT(IJK_FACET), BLOCKED_CELL_AT(IJK_FACET)
                     write(*,*)'PART POS',(DES_POS_NEW(LL,IBC),IBC=1,3)
                     write(*,*)'FACET NORM',(NORM_FACE(IBC,NF),IBC=1,3)
                     write(*,*)'BC_ID', BC_ID
                     write(*,*)'I,J,K (Facet)', I_FACET,J_FACET,K_FACET
                     call mfix_exit(1)
                  ENDIF
               ENDIF
               IF(FLUID_AT(IJK_FACET).AND.AREA>ZERO)THEN
                  DES_QW_Cond(IJK_FACET,phase_LL) = &
                     DES_QW_Cond(IJK_FACET, phase_LL) + QSWALL/AREA
               ENDIF

            ENDIF ! WALL BDRY
         ENDDO                  ! CELL_COUNT (facets)
      ENDDO  ! LL
      RETURN

   END SUBROUTINE CALC_DEM_THERMO_WITH_WALL_STL

END MODULE CALC_COLLISION_WALL
