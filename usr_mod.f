#include "error.inc"

MODULE usr
!
!     Declare the user-defined namelist variables (usrnlst.inc) in this module.
!     Also Include user-defined variables in this module.  To access the
!     variables from a subroutine add the statement "Use usr".  If allocatable
!     arrays are defined in this module allocate them in usr0.  To turn on the
!     user defined subroutines (usr0, usr1, and usr2) set keyword CALL_USR to true.

   use error_manager
   use compar, only: mype, pe_io
   use mpi_utility, only: bcast
   use DISCRETELEMENT
   use fs_util, only: file_exists
   use compar, only: mype, pe_io
   use exit, only: mfix_exit
   USE constant, only: pi
   USE param1, only: zero, one, half, undefined

   IMPLICIT NONE

! Maximum number of keyframes allowed
   integer, parameter :: maxcount = 100
! Keyframe flag
   logical, dimension(maxcount) ::  read_kf
! Total number of keyframe tables
   integer :: kf_count
! 插值速度
   DOUBLE PRECISION :: Current_V

! Struct for storing keyframe data
   type kf_struct
      integer :: ID
      integer :: nrows, nvars
      integer :: velocity_index
      character(len=512) :: filename
      character(len=10) :: interp_method
      double precision, allocatable :: velocity(:)
      double precision, allocatable :: var(:,:)
      double precision, allocatable :: var_current(:)
   end type kf_struct

   type(kf_struct), allocatable,  dimension(:) :: kf_data
   double precision, allocatable, dimension(:) :: interpolated_var
   integer, dimension(maxcount) :: list
   integer, dimension(maxcount) :: rlist

contains

subroutine change_LSD_parameters_particle(L,I)
!-----------------------------------------------
! Modules
!-----------------------------------------------
   USE discretelement
   IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! parameters of particle-particle contact pair    
   DOUBLE PRECISION, INTENT(IN) :: L,I
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------
! Particle-Particle Collision Parameters ------------------------------>
  
   ! Calculate masses used for collision calculations.
   MASS_L = REF_MASS(L)
   MASS_I = REF_MASS(I)
   MASS_EFF = MASS_I*MASS_L/(MASS_I+MASS_L)
   
   ! Calculate the M-L normal and tangential damping coefficients.
   IF(EN_Var(L) .NE. ZERO) THEN
      DES_ETAN(M,L) = 2.0D0*SQRT(KN*MASS_EFF) * ABS(LOG(EN_Var(L)))
      DES_ETAN(M,L) = DES_ETAN(M,L)/SQRT(PI*PI + (LOG(EN_Var(L))**2))
   ELSE
      DES_ETAN(M,L) = 2.0D0*SQRT(KN*MASS_EFF)
   ENDIF
   DES_ETAT(M,L) = DES_ETAT_FAC*DES_ETAN(M,L)
   
   ! Store the entries in the symmetric matrix.
   DES_ETAN(L,M) = DES_ETAN(M,L)
   DES_ETAT(L,M) = DES_ETAT(M,L)
   
   ! Calculate the collision time scale.
   TCOLL_TMP = PI/SQRT(KN/MASS_EFF - &
      ((DES_ETAN(M,L)/MASS_EFF)**2)/4.d0)
   TCOLL = MIN(TCOLL_TMP, TCOLL)

   RETURN
   
END SUBROUTINE change_LSD_parameters_particle

SUBROUTINE Cal_Rel_Vel(DES_VEL_L, DES_VEL_I, DES_RADIUS_L, DES_RADIUS_I, DES_OMEGA_L, DES_OMEGA_I, NORM, DIST_LI, V_REL_TRANS_NORM, VREL_T)

   !-----------------------------------------------
   ! Modules
   !-----------------------------------------------
         USE discretelement, only: DES_CROSSPRDCT
         IMPLICIT NONE
   !-----------------------------------------------
   ! Dummy arguments
   !-----------------------------------------------
   ! parameters of particle-particle contact pair    
         DOUBLE PRECISION, INTENT(IN) :: DES_VEL_L(3)
         DOUBLE PRECISION, INTENT(IN) :: DES_VEL_I(3)
         DOUBLE PRECISION, INTENT(IN) :: DES_OMEGA_L(3)
         DOUBLE PRECISION, INTENT(IN) :: DES_OMEGA_I(3)
         DOUBLE PRECISION, INTENT(IN) :: DES_RADIUS_L
         DOUBLE PRECISION, INTENT(IN) :: DES_RADIUS_I
   ! distance between particle centers
         DOUBLE PRECISION, INTENT(IN) :: DIST_LI
   ! unit normal vector along the line of contact pointing from
   ! particle L to particle II
         DOUBLE PRECISION, INTENT(IN) :: NORM(3)
   ! slip velocity at point of contact
         DOUBLE PRECISION, INTENT(OUT) :: VREL_T(3)
   ! normal component of relative contact velocity (scalar)
         DOUBLE PRECISION, INTENT(OUT) :: V_REL_TRANS_NORM
   !-----------------------------------------------
   ! Local variables
   !-----------------------------------------------
   ! translational relative velocity
         DOUBLE PRECISION :: VRELTRANS(3)
   ! rotational velocity at point of contact
         DOUBLE PRECISION :: V_ROT(3), OMEGA_SUM(3)
   ! distance from the contact point to the particle centers
         DOUBLE PRECISION :: DIST_CL, DIST_CI
   !-----------------------------------------------
   
   ! translational relative velocity
            VRELTRANS(:) = (DES_VEL_L - DES_VEL_I)
   
   ! calculate the distance from the particle center to the contact point,
   ! which is taken as the radical line
   ! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
            DIST_CL = (DIST_LI**2 + DES_RADIUS_L**2 - DES_RADIUS_I**2)/&
               (2.d0*DIST_LI)
            DIST_CI = DIST_LI - DIST_CL
   
            OMEGA_SUM(:) = DES_OMEGA_L*DIST_CL + DES_OMEGA_I*DIST_CI
   
   ! calculate the rotational relative velocity
         V_ROT = DES_CROSSPRDCT(OMEGA_SUM, NORM)
   
   ! total relative velocity
         VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)
   
   ! normal component of relative velocity (scalar)
         V_REL_TRANS_NORM = DOT_PRODUCT(VRELTRANS,NORM)
   
   ! slip velocity of the contact point
   ! Equation (8) in Tsuji et al. 1992
         VREL_T =  VRELTRANS(:) - V_REL_TRANS_NORM*NORM(:)
   
         RETURN
   
      END SUBROUTINE Cal_Rel_Vel


   SUBROUTINE Cal_Rel_Vel_Wall(DES_VEL, OMEGA,  NORM, DIST, VRN, VRT)
      IMPLICIT NONE
! Dummy arguments:
!---------------------------------------------------------------------//
! Magnitude of the total relative translational velocity.
      DOUBLE PRECISION, INTENT(OUT):: VRN
! Total relative translational velocity (vector).
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(OUT):: VRT
! Unit normal from particle center to closest point on stl (wall)
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: NORM, OMEGA, DES_VEL 
! Distance between particle center and stl (wall).
      DOUBLE PRECISION, INTENT(IN) :: DIST
! Local variables
!---------------------------------------------------------------------//
! Additional relative translational motion due to rotation
      DOUBLE PRECISION, DIMENSION(DIMN) :: V_ROT
! Total relative velocity at contact point
      DOUBLE PRECISION, DIMENSION(DIMN) :: VRELTRANS

! Total relative velocity + rotational contribution
         V_ROT = DIST * OMEGA
         V_ROT = DES_CROSSPRDCT(V_ROT, NORM)

      VRELTRANS(:) =  DES_VEL + V_ROT(:)

! magnitude of normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! total relative translational slip velocity at the contact point
! Equation (8) in Tsuji et al. 1992
      VRT(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
      END SUBROUTINE Cal_Rel_Vel_Wall


   subroutine change_LSD_parameters_Wall(M)
      IMPLICIT NONE
! Local Variables:
!---------------------------------------------------------------------//
      ! particle_global_ID
      integer :: M
      ! Collision length scale.
      DOUBLE PRECISION :: TCOLL, TCOLL_TMP
      ! Collision length scale.
      DOUBLE PRECISION :: MASS_M, MASS_L, MASS_EFF
      !......................................................................!

! Particle-Wall Collision Parameters ---------------------------------->

! Check particle-wall normal restitution coefficient.
      !  EN = DES_EN_WALL_INPUT(M)
! Calculate masses used for collision calculations.
      MASS_M = REF_MASS(M)
      MASS_EFF = MASS_M

! Calculate the M-Wall normal and tangential damping coefficients.
      IF(EN .NE. ZERO) THEN
         DES_ETAN_WALL(M) = 2.d0*SQRT(KN_W*MASS_EFF)*ABS(LOG(EN))
         DES_ETAN_WALL(M) = DES_ETAN_WALL(M)/SQRT(PI*PI+(LOG(EN))**2)
      ELSE
         DES_ETAN_WALL(M) = 2.D0*SQRT(KN_W*MASS_EFF)
      ENDIF
      DES_ETAT_WALL(M) = DES_ETAT_W_FAC*DES_ETAN_WALL(M)

! Calculate the collision time scale.
!       print *,DTSOLID
!       TCOLL_TMP = PI/SQRT(KN_W/MASS_EFF -                           &
!          ((DES_ETAN_WALL(M)/MASS_EFF)**2.d0)/4.d0)

! Store the smalled calculated collision time scale. This value is used
! in time-marching the DEM solids.
!       DTSOLID = TCOLL/DTSOLID_FACTOR
!       print *,DTSOLID
      
      end subroutine change_LSD_parameters


!----------------------------------------------------------------------!
! Below is the code for keyframe data management.                      !
! It is anticipated that users won't need to modify it.                !
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: read_keyframe_data                                      !
!  Author: Sathish Sanjeevi                         Date:23-Oct-2019   !
!                                                                      !
!  Purpose: Keyframes are data tables provided by user to perform user-!
!           defined operations such as ramping or other velocity-dependent !
!           operations. This functions reads keyframe files provided   !
!           by the user.                                               !
!                                                                      !
!----------------------------------------------------------------------!

   subroutine read_keyframe_data

! Dummy arguments:
!---------------------------------------------------------------------
      integer i, j, k
      integer id, fileunit, allocstatus, io
      character(len = 8) :: formt, xi
      double precision :: velocity_init


! Read and store keyframe data in root process
      if(myPE == pe_io) then


         i = 0
         list = 0
         rlist = 0

! Dry run to count actual number of keyframes
         do id = 1, maxcount
            if(read_kf(id) .eqv. .true.) then
               i = i+1

! Keyframe IDs can be non-contiguous. Efficient storage needs contiguous
! structs. Store indices based on KF IDs and vice versa.
               list(id) = i
               rlist(i) = id
            end if
         end do
         kf_count = i
         print *, "Total number of keyframe files : ", kf_count
         allocate(kf_data(kf_count), stat = allocstatus)

! Loop to read each keyframe file
         do i = 1, kf_count
            id = rlist(i)
            write(xi, '(I0.4)') id
            kf_data(i)%filename = 'data_kf_'//trim(xi)//'.txt'
            print *, 'Current keyframe file is : ', kf_data(i)%filename

! Report error if file does not exists
            if (.not. file_exists(kf_data(i)%filename)) then
               write(err_msg, 9999)  trim(kf_data(i)%filename)
9999           format(2/,1x,70('*')/' from: usr_mod.f',/' error 9999: ',    &
                  'input file does not exist: ',a/,'aborting.',/1x,   &
                  70('*'),2/)
               call log_error()
            endif

! Read keyframe file
            open(newunit = fileunit, file = kf_data(i)%filename)
            read(fileunit, *) kf_data(i)%nrows, kf_data(i)%nvars
            read(fileunit, *) kf_data(i)%interp_method
            read(fileunit, *)   ! Dummy read to ignore header
            kf_data(i)%velocity_index = 1
            allocate(kf_data(i)%velocity(kf_data(i)%nrows))
            allocate(kf_data(i)%var(kf_data(i)%nrows, &
               kf_data(i)%nvars))
! Actual reading of keyframe table
            do j = 1, kf_data(i)%nrows
               read(fileunit, *, iostat = io) kf_data(i)%velocity(j), &
                  (kf_data(i)%var(j,k), k = 1, kf_data(i)%nvars)
               ! check if file is in proper condition
               if (io>0) then
                  write(err_msg, 1111) j+3, kf_data(i)%filename
1111              format('Check if data is a valid floating point at line : ', i5/,&
                     'in file : ', a)
                  call log_error()
               else if (io<0) then
                  write(err_msg, 2222) kf_data(i)%filename
2222              format('End of file reached earlier in file : ', a)
                  call log_error()
               end if
            end do
            close(fileunit)


! Important: check if keyframe velocity are sorted in ascending order
            velocity_init = 0.D0
            do j = 1, kf_data(i)%nrows
               if(kf_data(i)%velocity(j) .ge. velocity_init) then
                  velocity_init = kf_data(i)%velocity(j)
               else
                  write(err_msg, 9555) trim(kf_data(i)%filename), j+2
9555              format('Velocity cannot be decreasing. Check keyframe file : ',a/,&
                     'first column (velocity) at line : ', i5)
                  call log_error()
               end if
            end do

! Initialize starting keyframe variable to the first one
            kf_data(i)%var_current = kf_data(i)%var(1,:)
            kf_data(i)%velocity_index = 1
         end do

         print *, "Keyframe read complete!"

      end if

! Broadcast information to all processes from root
      call bcast(kf_count, pe_io)
      call bcast(list, pe_io)
      call bcast(rlist, pe_io)
   end subroutine read_keyframe_data



!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: interpolate_keyframe_data                               !
!  Author: Sathish Sanjeevi                         Date:23-Oct-2019   !
!                                                                      !
!  Purpose: Interpolate keyframe data based on current simulation      !
!           velocitystep                                                   !
!                                                                      !
!----------------------------------------------------------------------!

   subroutine interpolate_keyframe_data(velocity)
      double precision, intent(in) :: velocity
      integer i
      integer vindex, nrows
      character(len=10) :: interp
      double precision :: vmin, vmax
      double precision :: kfvelocity, kfvelocity_next
      double precision :: denom
      double precision, allocatable :: kfvar(:), kfvar_next(:)
      logical :: flag

! Interpolate keyframe value based on current velocitystep
      do i = 1, kf_count
         kf_data(i)%velocity_index = 1
         nrows  = kf_data(i)%nrows
         interp = kf_data(i)%interp_method
         vmin = kf_data(i)%velocity(1)
         vmax   = kf_data(i)%velocity(nrows)
         vindex = kf_data(i)%velocity_index
         kfvelocity = kf_data(i)%velocity(vindex)
         kfvelocity_next = kf_data(i)%velocity(vindex+1)
         kfvar = kf_data(i)%var(vindex,:)
         kfvar_next = kf_data(i)%var(vindex+1,:)


! Interpolate only if the simulation velocity lies within bounds of
! keyframe velocity data
         if(velocity.gt.vmin.and.velocity.lt.vmax) then
            flag = velocity.gt.kf_data(i)%velocity(vindex).and.&
               velocity.le.kf_data(i)%velocity(vindex+1)

! Loop repeatedly until the simulation velocity falls within the
! corresponding keyframe velocity intervals. This can happen if keyframe
! velocitysteps are smaller than the simulation velocitysteps
            do while(.not.flag)
               vindex = vindex + 1
               kf_data(i)%velocity_index = vindex
               kfvelocity = kf_data(i)%velocity(vindex)
               kfvelocity_next = kf_data(i)%velocity(vindex+1)
               kfvar = kf_data(i)%var(vindex,:)
               kfvar_next = kf_data(i)%var(vindex+1,:)
               flag = velocity.gt.kf_data(i)%velocity(vindex).and.&
                  velocity.le.kf_data(i)%velocity(vindex+1)
            end do

! Interpolation can be a `linear` one or `step` (floor) function
            if(velocity.gt.kfvelocity.and.velocity.le.kfvelocity_next) then
               if(interp.eq.'linear') then
                  denom = kfvelocity_next - kfvelocity
! Check if denominator is zero
                  if(denom.gt.0.D0) then
                     kf_data(i)%var_current(:) = kfvar(:) + &
                        (velocity - kfvelocity)/denom * (kfvar_next(:) - kfvar(:))
                  else
                     kf_data(i)%var_current(:) = kfvar(:)
                  end if
               elseif(interp.eq.'step') then
                  kf_data(i)%var_current(:) = kfvar(:)
               else
                  write(*,5555)
5555              FORMAT('Enter a valid interpolation scheme: linear or step')
               end if
            end if
         end if

! Use the first or last keyframe velocity when simulation velocity is outside
! keyframe bounds
         if(velocity.le.vmin) &
            kf_data(i)%var_current(:) = kf_data(i)%var(1,:)
         if(velocity.ge.vmax) &
            kf_data(i)%var_current(:) = kf_data(i)%var(nrows,:)
      end do
   end subroutine interpolate_keyframe_data

END MODULE usr
