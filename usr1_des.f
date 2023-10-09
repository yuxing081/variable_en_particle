!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS1_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms have been calculated but before they are     !
!  applied. The user may insert code in this routine or call user      !
!  defined subroutines.                                                !
!                                                                      !
!  This routine is called from the time loop, but no indices (fluid    !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR1_DES(Current_V,L,I)

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! parameters of particle-particle contact pair    
      DOUBLE PRECISION, INTENT(IN) :: L,I,Current_V
      integer :: M

      M = 1

      call interpolate_keyframe_data(Current_V)
      EN_Var(L) = kf_data(M)%var_current(1)
      call change_LSD_parameters_particle(L,I)
      
      RETURN
      END SUBROUTINE USR1_DES
