!+---------------------------------------------------------------------+
!| this module contains variables, arraies and subroutines related to  |
!| the mhd functionality.                                              |
!+---------------------------------------------------------------------+
!| change record                                                       |
!| -------------                                                       |
!| 28-Oct-2022  | Created by J. Fang STFC Daresbury Laboratory         |
!+---------------------------------------------------------------------+
module mhd
  !
  use decomp_2d, only : mytype,xsize,nrank,nproc
  use hdf5
  use h5lt
  use partack, only: mpistop,rankname
  !
  implicit none
  !
  logical :: mhd_active
  real(8) :: hartmann,stuart 
  !+------------+------------------------------------------------------+
  !|  mhd_active| the swith to activate the mhd module.                |
  !|    hartmann| hartmann number, the ratio of lorentz force to       |
  !|            | viscous force                                        |
  !|      stuart| Stuart number, magnetic interaction parameter, ratio |
  !|            | of electromagnetic to inertial forces                |
  !+------------+------------------------------------------------------+
  !
  real(mytype),allocatable,dimension(:,:,:,:) :: magent,magelf,elefld
  real(mytype),allocatable,dimension(:,:,:) :: elcpot
  !+------------+------------------------------------------------------+
  !|      magent| magnetic field                                       |
  !|      magelf| Electromagnetic forced to apply to the momentum eq.  |
  !|      elcpot| electric potential                                   |
  !|      elefld| electric field as the gradient of potential          |
  !+------------+------------------------------------------------------+
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is to initilise the MHD module.                   |
  !+-------------------------------------------------------------------+
  !| change record                                                     |
  !| -------------                                                     |
  !| 28-Oct-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine mhd_init
    !
    use param, only: re
    !
    stuart=hartmann**2/re
    !
    if(nrank==0) then
      !
      if(mhd_active) then
        print*,'** MHD Module activated'
        print*,'** Hartmann number: ',hartmann
        print*,'** Reynolds number: ',re
        print*,'** Stuart number  : ',stuart
      endif
      !
    endif
    !
    allocate( magent(xsize(1),xsize(2),xsize(3),1:3),                  &
              magelf(xsize(1),xsize(2),xsize(3),1:3),                  &
              elefld(xsize(1),xsize(2),xsize(3),1:3),                  &
              elcpot(xsize(1),xsize(2),xsize(3)) )
    !
    if(nrank==0) print*,'** MHD fields allocated'
    !
    magent(:,:,:,1)=0.d0
    magent(:,:,:,2)=1.d0
    magent(:,:,:,3)=0.d0
    !
    if(nrank==0) print*,'** magnetic field initilised'
    !
  end subroutine mhd_init
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mhd_init.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is used to add the electromagnetic force to the   |
  !| momentum equation                                                 |
  !+-------------------------------------------------------------------+
  !| change record                                                     |
  !| -------------                                                     |
  !| 28-Oct-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine momentum_forcing_mhd(dux1,duy1,duz1,ux1,uy1,uz1)
    !
    use param,     only : dx,dz
    use variables, only : yp,ny,nz
    use decomp_2d, only : xstart
    use constants, only : pi
    !
    ! arguments
    real(mytype),intent(in), dimension(xsize(1),xsize(2),xsize(3)) ::  &
                                                ux1,uy1,uz1
    real(mytype),intent(inout),                                        &
                dimension(xsize(1),xsize(2),xsize(3)) ::  dux1,duy1,duz1
    !
    ! local data
    integer :: i,j,k
    real(mytype) :: elecur(3),var1(3),var2(3)
    !
    ! real(mytype) :: xx(xsize(1)),yy(xsize(2)),zz(xsize(3))
    ! !
    ! do i=1,xsize(1)
    !   xx(i)=real(i-1,mytype)*dx
    ! enddo
    ! do j=1,xsize(2)
    !   !
    !   if(j+xstart(2)-1>ny) then
    !     yy(j)=2.d0*yp(ny)-yp(ny-1)
    !   elseif(j+xstart(2)-1<1) then
    !     yy(j)=2.d0*yp(1)-yp(2)
    !   else
    !     yy(j)=yp(j+xstart(2)-1)
    !   endif
    !   !
    ! enddo
    ! do k=1,xsize(3)
    !   zz(k)=real((k+xstart(3)-2),mytype)*dz
    ! enddo
    ! !
    ! do k = 1, xsize(3)
    ! do j = 1, xsize(2)
    ! do i = 1, xsize(1)
    !   ! elcpot(i,j,k)=sin(xx(i)*pi)
    !   ! elcpot(i,j,k)=sin(zz(k)*pi)
    !   elcpot(i,j,k)=sin(0.5d0*yy(j)*pi)
    ! enddo
    ! enddo
    ! enddo
    !
    elefld=electric_field(elcpot)
    ! to obtain the electric field from electric potential
    !
    ! write(rankname,'(i4.4)')nrank
    ! open(18,file='testout/profile'//rankname//'.dat')
    ! ! do i = 1, xsize(1)
    !   ! write(18,*)xx(i),elcpot(i,1,1),elefld(i,1,1,1)
    ! ! do k = 1, xsize(3)
    ! !   write(18,*)zz(k),elcpot(1,1,k),elefld(1,1,k,3)
    ! do j = 1, xsize(2)
    !   write(18,*)yy(j),elcpot(1,j,1),elefld(1,j,1,2)
    ! enddo
    ! close(18)
    ! !
    !
    do k = 1, xsize(3)
    do j = 1, xsize(2)
    do i = 1, xsize(1)
      !
      var1(1)=ux1(i,j,k)
      var1(2)=uy1(i,j,k)
      var1(3)=uz1(i,j,k)
      !
      var2=cross_product(var1,magent(i,j,k,:))
      !
      elecur(:)=elefld(i,j,k,:)+var2
      !
      magelf(i,j,k,:)=cross_product(elecur,magent(i,j,k,:))*stuart
      !
      dux1(i,j,k) = dux1(i,j,k)+magelf(i,j,k,1)
      duy1(i,j,k) = duy1(i,j,k)+magelf(i,j,k,2)
      duz1(i,j,k) = duz1(i,j,k)+magelf(i,j,k,3)
      !
    enddo
    enddo
    enddo
    !
    call mpistop
    !
  end subroutine momentum_forcing_mhd
  !+-------------------------------------------------------------------+
  !| The end of the subroutine momentum_forcing.                       |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is used to calculate the electric field from      |
  !| electric potential.                                               |
  !+-------------------------------------------------------------------+
  !| change record                                                     |
  !| -------------                                                     |
  !| 28-Oct-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  function electric_field(epotential) result(efield)
    !
    use decomp_2d, only : ysize,zsize,transpose_x_to_y,                &
                          transpose_y_to_z,transpose_y_to_x,           &
                          transpose_z_to_y
    use variables, only: ffxpS,fsxpS,fwxpS,ffypS,fsypS,fwypS,ffzpS,    &
                         fszpS,fwzpS,ppy,sx,sy,sz,derxs,derys,derzs
    use param, only: zero
    use var, only : ta1,di1,ta2,di2,ta3,di3,td1,td2,td3,tg1
    !
    real(8) :: efield(xsize(1),xsize(2),xsize(3),3)
    real(8),intent(in) :: epotential(xsize(1),xsize(2),xsize(3))
    !
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: dpot1
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: pot2, dpot2
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: pot3, dpot3
    !
    call transpose_x_to_y(epotential,pot2)
    call transpose_y_to_z(      pot2,pot3)
    !
    call derxS (dpot1, epotential, di1, sx, ffxpS, fsxpS, fwxpS,      xsize(1), xsize(2), xsize(3), 1, zero)
    call deryS (dpot2,       pot2, di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
    call derzS (dpot3,       pot3, di3, sz, ffzpS, fszpS, fwzpS,      zsize(1), zsize(2), zsize(3), 1, zero)
    !
    efield(:,:,:,1)=-dpot1
    !
    call transpose_y_to_x(dpot2,dpot1)
    !
    efield(:,:,:,2)=-dpot1
    !
    call transpose_z_to_y(dpot3,dpot2)
    call transpose_y_to_x(dpot2,dpot1)
    
    efield(:,:,:,3)=-dpot1
    !
    return
    !
  end function electric_field
  !+-------------------------------------------------------------------+
  !| The end of the subroutine electric_field.                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is used to solve the poisson equation related to  |
  !| MHD.                                                              |
  !+-------------------------------------------------------------------+
  !| change record                                                     |
  !| -------------                                                     |
  !| 28-Oct-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine solve_mhd_poisson(ux1,uy1,uz1)

    use decomp_2d, only : mytype, xsize, zsize, ph1, nrank
    use decomp_2d_poisson, only : poisson
    use var, only : nzmsize
    use var, only : dv3
    use param, only : ntime, nrhotime, npress
    use param, only : ilmn, ivarcoeff, zero, one 

    implicit none

    !! inputs
    real(mytype),dimension(xsize(1), xsize(2), xsize(3)),intent(in) :: ux1, uy1, uz1
    !
  end subroutine solve_mhd_poisson
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solve_mhd_poisson.                      |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to do the cross product for a 3-D vector.        | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  pure function cross_product(a,b)
    !
    real(8) :: cross_product(3)
    real(8),dimension(3),intent(in) :: a(3), b(3)
  
    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)
    !
    return
    !
  end function cross_product
  !+-------------------------------------------------------------------+
  !| The end of the function cross_product.                            |
  !+-------------------------------------------------------------------+
  !
end module mhd
!+---------------------------------------------------------------------+
! the end of the module mhd                                            |
!+---------------------------------------------------------------------+