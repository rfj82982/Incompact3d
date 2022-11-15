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
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) ::   &
                                                ux1,uy1,uz1
    real(mytype),intent(inout),                                        &
                dimension(xsize(1),xsize(2),xsize(3)) ::  dux1,duy1,duz1
    !
    real(mytype) :: eforce(3)
    ! local data
    integer :: i,j,k
    real(mytype) :: elecur(3),var1(3),var2(3)
    !
    real(mytype) :: xx(xsize(1)),yy(xsize(2)),zz(xsize(3))
    !
    !
    call solve_mhd_poisson(ux1,uy1,uz1)
    !
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
    ! elefld=electric_field(elcpot)
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
    ! write(rankname,'(i4.4)')nrank
    ! open(18,file='mhd_force'//rankname//'.dat')
    !
    do k = 1, xsize(3)
    do j = 1, xsize(2)
    do i = 1, xsize(1)
      !
      eforce=cross_product(elefld(i,j,k,:),magent(i,j,k,:))*stuart
      !
      dux1(i,j,k) = dux1(i,j,k)+eforce(1)
      duy1(i,j,k) = duy1(i,j,k)+eforce(2)
      duz1(i,j,k) = duz1(i,j,k)+eforce(3)
      !
      ! if(i==1 .and. k==1) then
      !   write(18,*)yy(j),magent(i,j,k,2),elefld(i,j,k,3),eforce(1)
      ! endif
      !
    enddo
    enddo
    enddo
    !
    ! close(18)
    ! print*,' << mhd_force',rankname,'.dat'
    ! !
    ! call mpistop
    !
  end subroutine momentum_forcing_mhd
  !+-------------------------------------------------------------------+
  !| The end of the subroutine momentum_forcing.                       |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is used to calculate the gradient of a general    |
  !| function on velocity mesh..................                       |
  !+-------------------------------------------------------------------+
  !| change record                                                     |
  !| -------------                                                     |
  !| 28-Oct-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  function grad_vmesh(phi) result(dphi)
    !
    use decomp_2d, only : ysize,zsize,transpose_x_to_y,                &
                          transpose_y_to_z,transpose_y_to_x,           &
                          transpose_z_to_y
    use variables, only: ffxpS,fsxpS,fwxpS,ffypS,fsypS,fwypS,ffzpS,    &
                         fszpS,fwzpS,ppy,sx,sy,sz,derxs,derys,derzs
    use param, only: zero
    use var, only : ta1,di1,ta2,di2,ta3,di3,td1,td2,td3,tg1
    !
    real(mytype) :: dphi(xsize(1),xsize(2),xsize(3),3)
    real(mytype),intent(in) :: phi(xsize(1),xsize(2),xsize(3))
    !
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: dpot1
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: pot2, dpot2
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: pot3, dpot3
    !
    call transpose_x_to_y(phi,pot2)
    call transpose_y_to_z(      pot2,pot3)
    !
    call derxS (dpot1,  phi, di1, sx, ffxpS, fsxpS, fwxpS,      xsize(1), xsize(2), xsize(3), 1, zero)
    call deryS (dpot2, pot2, di2, sy, ffypS, fsypS, fwypS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)
    call derzS (dpot3, pot3, di3, sz, ffzpS, fszpS, fwzpS,      zsize(1), zsize(2), zsize(3), 1, zero)
    !
    dphi(:,:,:,1)=dpot1
    !
    call transpose_y_to_x(dpot2,dpot1)
    !
    dphi(:,:,:,2)=dpot1
    !
    call transpose_z_to_y(dpot3,dpot2)
    call transpose_y_to_x(dpot2,dpot1)
    
    dphi(:,:,:,3)=dpot1
    !
    return
    !
  end function grad_vmesh
  !+-------------------------------------------------------------------+
  !| The end of the subroutine grad_vmesh.                             |
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
    use var, only : nzmsize,dv3
    use param, only : ntime, nrhotime, npress,ilmn, ivarcoeff, zero, one 
    use navier,only : gradp

    implicit none

    !! inputs
    real(mytype),dimension(xsize(1), xsize(2), xsize(3)),intent(in) :: ux1, uy1, uz1
    !
    !! local data
    real(mytype),dimension(xsize(1), xsize(2), xsize(3), 3) :: ucB
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: rhs
    
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: div3
    !
    integer :: i,j,k,nlock
    real(mytype) :: var1(3),var2(3)
    logical :: converged
    !
    nlock=1
    !
    do k = 1, xsize(3)
    do j = 1, xsize(2)
    do i = 1, xsize(1)
      !
      var1(1)=ux1(i,j,k)
      var1(2)=uy1(i,j,k)
      var1(3)=uz1(i,j,k)
      !
      ucB(i,j,k,:) =cross_product(var1,magent(i,j,k,:))
      !
    enddo
    enddo
    enddo
    !
    rhs=divergence_sclar(ucB,nlock)
    !
    converged=.false.
    !
    do while(.not.converged)
      !
      call poisson(rhs)
      !
      CALL gradp(elefld(:,:,:,1),elefld(:,:,:,2),elefld(:,:,:,3),rhs)
      !
      converged=.true.
    enddo
    !
    do k = 1, xsize(3)
    do j = 1, xsize(2)
    do i = 1, xsize(1)
      !
      elefld(i,j,k,:)=-elefld(i,j,k,:)+ucB(i,j,k,:)
      !
    enddo
    enddo
    enddo
    !
    ! div3=divergence_sclar(elefld,nlock)
    ! !
    ! call div_check(div3)
    !
    ! do k = ph1%zst(1),ph1%zen(1)
    ! do j = ph1%zst(2),ph1%zen(2)
    ! do i = 1,nzmsize
    !   if(nrank==0) then
    !     print*,i,j,k,div3(i,j,k)
    !   endif
    ! enddo
    ! enddo
    ! enddo
    ! print*,xsize(1), xsize(2), xsize(3)
    ! print*,ph1%zst(1),ph1%zen(1),ph1%zst(2),ph1%zen(2),nzmsize
    ! !
    ! call mpistop
    !
  end subroutine solve_mhd_poisson
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solve_mhd_poisson.                      |
  !+-------------------------------------------------------------------+
  !
  subroutine  div_check(divec)
    !
    USE decomp_2d
    use partack, only: psum,pmax
    !
    real(mytype),intent(in) :: divec(:,:,:)
    !
    integer :: i,j,k
    real(mytype) :: tmax,tmoy
    !
    tmax=-1609._mytype
    tmoy=0._mytype
    do k = 1, size(divec,3)
    do j = 1, size(divec,2)
    do i = 1, size(divec,1)
      if (divec(i,j,k).gt.tmax) tmax=divec(i,j,k)
      tmoy=tmoy+abs(divec(i,j,k))
    enddo
    enddo
    enddo
    tmoy=tmoy/(size(divec,1)*size(divec,2)*size(divec,3))
    !
    ! call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
    ! call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    tmax=pmax(tmax)
    tmoy=psum(tmoy)
    !
    if (nrank == 0) then
       write(*,*) 'DIV J max mean=',real(tmax,mytype),real(tmoy/real(nproc),mytype)
    endif
    !
  end subroutine
  !
  subroutine  divergence_vmesh_check(vec,divec)
    !
    USE MPI
    USE decomp_2d
    use partack, only: psum,pmax
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),3),intent(in) :: vec
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(out) :: divec
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),3) :: dvec1,dvec2,dvec3
    real(mytype) :: div,tmax,tmoy,tmax1,tmoy1,code
    !
    integer :: i,j,k
    !
    dvec1=grad_vmesh(vec(:,:,:,1))
    dvec2=grad_vmesh(vec(:,:,:,2))
    dvec3=grad_vmesh(vec(:,:,:,3))
    !
    tmax=-1609._mytype
    tmoy=0._mytype
    do k = 1, xsize(3)
    do j = 1, xsize(2)
    do i = 1, xsize(1)
      divec(i,j,k)=div
      !
      div=dvec1(i,j,k,1)+dvec2(i,j,k,2)+dvec3(i,j,k,3)
      if (div.gt.tmax) tmax=div
      tmoy=tmoy+abs(div)
    enddo
    enddo
    enddo
    tmoy=tmoy/(xsize(1)*xsize(2)*xsize(3))
    !
    ! call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
    ! call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    tmax=pmax(tmax)
    tmoy=psum(tmoy)
    !
    if (nrank == 0) then
       write(*,*) 'DIV J max mean=',real(tmax,mytype),real(tmoy/real(nproc),mytype)
    endif
    !
  end subroutine  divergence_vmesh_check
  !+-------------------------------------------------------------------+
  !| The end of the subroutine divergence_vmesh_check.                 |
  !+-------------------------------------------------------------------+
  !
  !!############################################################################
  !subroutine DIVERGENCe
  !Calculation of div u* for nlock=1 and of div u^{n+1} for nlock=2
  ! input :  vec (on velocity mesh)
  ! output : pp3 (on pressure mesh)
  !written by SL 2018
  !############################################################################
  function divergence_sclar(vec,nlock) result(pp3)

    USE param
    USE decomp_2d
    USE variables
    USE var, ONLY: ta1, tb1, tc1, pp1, pgy1, pgz1, di1, &
         duxdxp2, uyp2, uzp2, duydypi2, upi2, ta2, dipp2, &
         duxydxyp3, uzp3, po3, dipp3, nxmsize, nymsize, nzmsize
    USE MPI
    USE ibm_param
    use navier, only: extrapol_drhodt

    implicit none

    !  TYPE(DECOMP_INFO) :: ph1,ph3,ph4

    !X PENCILS NX NY NZ  -->NXM NY NZ
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),3),intent(in) :: vec
    !
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: pp3

    integer :: nvect3,i,j,k,nlock
    integer :: code
    real(mytype) :: tmax,tmoy,tmax1,tmoy1

    nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

    ta1(:,:,:) = vec(:,:,:,1)
    tb1(:,:,:) = vec(:,:,:,2)
    tc1(:,:,:) = vec(:,:,:,3)

    !WORK X-PENCILS

    call derxvp(pp1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)

    if (ilmn.and.(nlock.gt.0)) then
       call interxvp(pgy1,ta1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
       pp1(:,:,:) = pp1(:,:,:) + pgy1(:,:,:)
    endif

    call interxvp(pgy1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
    call interxvp(pgz1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

    call transpose_x_to_y(pp1,duxdxp2,ph4)!->NXM NY NZ
    call transpose_x_to_y(pgy1,uyp2,ph4)
    call transpose_x_to_y(pgz1,uzp2,ph4)

    !WORK Y-PENCILS
    call interyvp(upi2,duxdxp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
    call deryvp(duydypi2,uyp2,dipp2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)

    !! Compute sum dudx + dvdy
    duydypi2(:,:,:) = duydypi2(:,:,:) + upi2(:,:,:)

    call interyvp(upi2,uzp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

    call transpose_y_to_z(duydypi2,duxydxyp3,ph3)!->NXM NYM NZ
    call transpose_y_to_z(upi2,uzp3,ph3)

    !WORK Z-PENCILS
    call interzvp(pp3,duxydxyp3,dipp3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
    call derzvp(po3,uzp3,dipp3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,0)

    !! Compute sum dudx + dvdy + dwdz
    pp3(:,:,:) = pp3(:,:,:) + po3(:,:,:)

    if (nlock==2) then
       pp3(:,:,:)=pp3(:,:,:)-pp3(ph1%zst(1),ph1%zst(2),nzmsize)
    endif

    ! tmax=-1609._mytype
    ! tmoy=zero
    ! do k=1,nzmsize
    !    do j=ph1%zst(2),ph1%zen(2)
    !       do i=ph1%zst(1),ph1%zen(1)
    !          if (pp3(i,j,k).gt.tmax) tmax=pp3(i,j,k)
    !          tmoy=tmoy+abs(pp3(i,j,k))
    !       enddo
    !    enddo
    ! enddo
    ! tmoy=tmoy/nvect3

    ! call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
    ! call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    ! if ((nrank == 0) .and. (nlock > 0).and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
    !    if (nlock == 2) then
    !       write(*,*) 'DIV ucB  max mean=',real(tmax1,mytype),real(tmoy1/real(nproc),mytype)
    !    else
    !       write(*,*) 'DIV ucB* max mean=',real(tmax1,mytype),real(tmoy1/real(nproc),mytype)
    !    endif
    ! endif

    return
    !
  end function divergence_sclar
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