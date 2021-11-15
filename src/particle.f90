!+---------------------------------------------------------------------+
!| This module contains variables, arraies and subroutines related to  |
!| particle tracking.                                                  |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 15-11-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module particle
  !
  use decomp_2d, only : mytype
  !
  implicit none
  !
  ! particles
  logical :: lpartack
  integer :: numparticle
  real(mytype),allocatable,dimension(:) :: xpa,ypa,zpa
  real(mytype),allocatable,dimension(:) :: ux_pa,uy_pa,uz_pa
  !+------------------+--------------------------------------------+
  !|         lpartack | switch of particel tracking                |
  !|      numparticle | number of particles in the domain          |
  !|      xpa,ypa,zpa | x,y,z coordinates of particles             |
  !|            ux_pa |                                            |
  !|            uy_pa |                                            |
  !|            uz_pa | velocity of particles                      |
  !+------------------+--------------------------------------------+
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to allocate particle arraies.                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine partialloc
    !
    !! particle tracking
    allocate(xpa(1:numparticle),ypa(1:numparticle),zpa(1:numparticle))
    allocate(ux_pa(1:numparticle),uy_pa(1:numparticle),uz_pa(1:numparticle))
    !
  end subroutine partialloc
  !+-------------------------------------------------------------------+
  ! The end of the subroutine partialloc                               |
  !+-------------------------------------------------------------------+
    !
end module particle
!+---------------------------------------------------------------------+
! The end of the module particle                                       |
!+---------------------------------------------------------------------+