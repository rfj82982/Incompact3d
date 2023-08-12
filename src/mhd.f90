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
  use mptool, only: mpistop,rankname,psum,pmax,cross_product
  !
  implicit none
  !
  logical :: mhd_active,mhd_equation,sync_Bm_needed= .true.
  real(8) :: hartmann,stuart,rem
  !+------------+------------------------------------------------------+
  !|  mhd_active| the swith to activate the mhd module.                |
  !|    hartmann| hartmann number, the ratio of lorentz force to       |
  !|            | viscous force                                        |
  !|      stuart| Stuart number, magnetic interaction parameter, ratio |
  !|            | of electromagnetic to inertial forces                |
  !+------------+------------------------------------------------------+
  !
  real(mytype),allocatable,dimension(:,:,:,:) :: Bm,magelf,Je
  real(mytype),allocatable,dimension(:,:,:,:,:) :: dBm
  real(mytype),allocatable,dimension(:,:,:) :: elcpot
  !+------------+------------------------------------------------------+
  !|          Bm| magnetic field                                       |
  !|      magelf| Electromagnetic forced to apply to the momentum eq.  |
  !|      elcpot| electric potential                                   |
  !|          Je| electric field as the gradient of potential          |
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
    use param, only: re,ntime
    !
    ! stuart=hartmann**2/re
    !
    if(stuart<=1.d-15) then
      stuart=hartmann**2/re
    endif
    if(hartmann<=1.d-15) then
      hartmann=sqrt(stuart*re)
    endif
    !
    if(nrank==0) then
      !
      print*,'** MHD Module activated'
      print*,'**    MHD equation:',mhd_equation
      print*,'** Hartmann number: ',hartmann
      print*,'**   Stuart number: ',stuart
      print*,'**              Re: ',re
      print*,'**     Magnetic Re: ',rem
      !
    endif
    !
    allocate( Bm(xsize(1),xsize(2),xsize(3),1:3),                  &
              magelf(xsize(1),xsize(2),xsize(3),1:3),              &
              Je(xsize(1),xsize(2),xsize(3),1:3),                  &
              elcpot(xsize(1),xsize(2),xsize(3)),                  &
              dBm(xsize(1),xsize(2),xsize(3),1:3,1:ntime) )
    !
    if(nrank==0) print*,'** MHD fields allocated'
    !
    ! Bm(:,:,:,1)=0.d0
    ! Bm(:,:,:,2)=1.d0
    ! Bm(:,:,:,3)=0.d0
    ! !
    ! if(nrank==0) print*,'** magnetic field initilised'
    ! !
  end subroutine mhd_init
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mhd_init.                               |
  !+-------------------------------------------------------------------+
  !
  subroutine int_time_magnet
    !
    USE param
    USE variables
    USE decomp_2d
    use ydiff_implicit

    implicit none
    !
    integer :: k

    if( iimplicit == 1 ) then
        call inttimp(Bm(:,:,:,1), dBm(:,:,:,1,:), 0, -1, mhdvar=1 )
        call inttimp(Bm(:,:,:,2), dBm(:,:,:,2,:), 0, -1, mhdvar=2 )
        call inttimp(Bm(:,:,:,3), dBm(:,:,:,3,:), 0, -1, mhdvar=3 )
    else
       if(itimescheme.eq.3) then
           !>>> Adams-Bashforth third order (AB3)

           ! Do first time step with Euler
           if(itime.eq.1.and.irestart.eq.0) then
              Bm=dt*dBm(:,:,:,:,1)+Bm
           elseif(itime.eq.2.and.irestart.eq.0) then
              ! Do second time step with AB2
              Bm=onepfive*dt*dBm(:,:,:,:,1)-half*dt*dBm(:,:,:,:,2)+Bm
              dBm(:,:,:,:,3)=dBm(:,:,:,:,2)
           else
              ! Finally using AB3
              Bm=adt(itr)*dBm(:,:,:,:,1)+bdt(itr)*dBm(:,:,:,:,2)+cdt(itr)*dBm(:,:,:,:,3)+Bm
              dBm(:,:,:,:,3)=dBm(:,:,:,:,2)
           endif
           dBm(:,:,:,:,2)=dBm(:,:,:,:,1)

        elseif(itimescheme.eq.5) then
          !
           if(itr.eq.1) then
              Bm=gdt(itr)*dBm(:,:,:,:,1)+Bm
           else
              Bm=adt(itr)*dBm(:,:,:,:,1)+bdt(itr)*dBm(:,:,:,:,2)+Bm
           endif
           dBm(:,:,:,:,2)=dBm(:,:,:,:,1)
           !
       endif
   endif
    !
  end subroutine int_time_magnet
  !
  function vortcal(dux1,duy1,duz1) result(omega)
    !
    USE decomp_2d
    USE variables
    USE param
    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use ibm_param, only : ubcx,ubcy,ubcz
    !
    ! real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: omega
    !
    ! Perform communications if needed
    ! if (sync_vel_needed) then
    !   call transpose_x_to_y(ux1,ux2)
    !   call transpose_x_to_y(uy1,uy2)
    !   call transpose_x_to_y(uz1,uz2)
    !   call transpose_y_to_z(ux2,ux3)
    !   call transpose_y_to_z(uy2,uy3)
    !   call transpose_y_to_z(uz2,uz3)
    !   sync_vel_needed = .false.
    ! endif

    ! !x-derivatives
    ! call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    ! call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    ! call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    ! !y-derivatives
    ! call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    ! call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    ! call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    ! !!z-derivatives
    ! call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    ! call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    ! call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    ! !!all back to x-pencils
    ! call transpose_z_to_y(ta3,td2)
    ! call transpose_z_to_y(tb3,te2)
    ! call transpose_z_to_y(tc3,tf2)
    ! call transpose_y_to_x(td2,tg1)
    ! call transpose_y_to_x(te2,th1)
    ! call transpose_y_to_x(tf2,ti1)
    ! call transpose_y_to_x(ta2,td1)
    ! call transpose_y_to_x(tb2,te1)
    ! call transpose_y_to_x(tc2,tf1)

    ! omega(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
    !                   + (tg1(:,:,:)-tc1(:,:,:))**2 &
    !                   + (tb1(:,:,:)-td1(:,:,:))**2)
    omega(:,:,:)=sqrt(  (duz1(:,:,:,2)-duy1(:,:,:,3))**2 &
                      + (dux1(:,:,:,3)-duz1(:,:,:,1))**2 &
                      + (duy1(:,:,:,1)-dux1(:,:,:,2))**2)
    !
    return
    !
  end function vortcal
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is calculate and output statistics of MHD flow.   |
  !+-------------------------------------------------------------------+
  !| change record                                                     |
  !| -------------                                                     |
  !| 01-May-2023  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine mhd_sta(ux1,uy1,uz1,dux1,duy1,duz1)
    !
    use decomp_2d
    use param,     only : ntime,t,nclx1, ncly1, nclz1,re
    use var,       only : itime
    use variables, only : nx, ny, nz, nxm, nym, nzm
    use mptool,    only : pmax,psum
    !
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    !
    ! local data
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: vorticity
    real(mytype) :: Ek,Em,Omegak,Omegam,Omgmax,Jmax,Omega(3),var1,var2,disrat
    logical,save :: lfirstcal=.true.
    integer,save :: nxc,nyc,nzc
    integer :: i,j,k
    !
    if(lfirstcal) then
      !
      if(nrank==0) then
        open(13,file='flowstat.dat')
        write(13,"(A7,1X,A13,7(1X,A20))")'itime','time',              &
                                'Ek','Em','Ωk','Ωm','ε','Ωmax','Jmax'

      endif
      !
      if (nclx1==1) then
         nxc=nxm
      else
         nxc=nx
      endif
      if (ncly1==1) then
         nyc=nym
      else
         nyc=ny
      endif
      if (nclz1==1) then
         nzc=nzm
      else
         nzc=nz
      endif
      !
      lfirstcal=.false.
      !
    endif
    !
    vorticity=vortcal(dux1,duy1,duz1)
    !
    Ek=0.d0
    Em=0.d0
    Omegak=0.d0
    Omegam=0.d0
    Omgmax=0.d0
    Jmax=0.d0
    do k=1,xsize(3)
    do j=1,xsize(2)
    do i=1,xsize(1)

      var1=vorticity(i,j,k)**2 
      var2=Je(i,j,k,1)**2+Je(i,j,k,2)**2+Je(i,j,k,3)**2

      Ek    =Ek    + ux1(i,j,k)**2+uy1(i,j,k)**2+uz1(i,j,k)**2
      Em    =Em    + Bm(i,j,k,1)**2+Bm(i,j,k,2)**2+Bm(i,j,k,3)**2
      Omegak=Omegak+ var1
      Omegam=Omegam+ var2
      Omgmax= max(Omgmax,var1)
      Jmax  = max(Jmax,var2)

    enddo
    enddo
    enddo
    !
    Ek    =psum(Ek    )
    Em    =psum(Em    )
    Omegak=psum(Omegak)
    Omegam=psum(Omegam)
    Omgmax=pmax(Omgmax)
    Jmax  =pmax(Jmax)
    !
    Ek    =Ek    /dble(nxc*nyc*nzc)/2.d0
    Em    =Em    /dble(nxc*nyc*nzc)/2.d0
    Omegak=Omegak/dble(nxc*nyc*nzc)/2.d0
    Omegam=Omegam/dble(nxc*nyc*nzc)/2.d0*Rem*Rem
    Omgmax=sqrt(Omgmax)
    Jmax  =sqrt(Jmax)*Rem
    !
    disrat=Ek/re+Em/rem
    ! print*,nxc,nyc,nzc
    !
    if(nrank==0) then
      write(13,"(I7,1X,E13.6E2,7(1X,E20.13E2))")itime,t,Ek,Em,Omegak, &
                                            Omegam,disrat,Omgmax,Jmax
    endif
    !
  end subroutine mhd_sta
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mhd_sta.                                |
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
    USE decomp_2d
    use mpi
    use param,     only : dx,dz
    use variables, only : ppy
    use variables, only : yp,ny,nz
    use decomp_2d, only : xstart
    use constants, only : pi
    use param, only : zero,two,three,dy,yly
    !
    ! arguments
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) ::   &
                                                ux1,uy1,uz1
    real(mytype),intent(inout),                                        &
                dimension(xsize(1),xsize(2),xsize(3)) ::  dux1,duy1,duz1
    !
    real(mytype) :: eforce(3),Ebar(3)
    real(mytype) :: ub,uball,coeff
    ! local data
    integer :: i,j,k,jloc,code
    real(mytype) :: elecur(3),var1(3),var2(3)
    !
    real(mytype) :: xx(xsize(1)),yy(xsize(2)),zz(xsize(3))
    !
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
    ! Je=electric_field(elcpot)
    ! to obtain the electric field from electric potential
    !
    ! write(rankname,'(i4.4)')nrank
    ! open(18,file='testout/profile'//rankname//'.dat')
    ! ! do i = 1, xsize(1)
    !   ! write(18,*)xx(i),elcpot(i,1,1),Je(i,1,1,1)
    ! ! do k = 1, xsize(3)
    ! !   write(18,*)zz(k),elcpot(1,1,k),Je(1,1,k,3)
    ! do j = 1, xsize(2)
    !   write(18,*)yy(j),elcpot(1,j,1),Je(1,j,1,2)
    ! enddo
    ! close(18)
    ! !
    ! write(rankname,'(i4.4)')nrank
    ! open(18,file='mhd_force'//rankname//'.dat')
    !
    ! ub = zero
    ! uball = zero
    ! coeff = dy / (yly * real(xsize(1) * zsize(3), kind=mytype))

    ! do k = 1, xsize(3)
    !    do jloc = 1, xsize(2)
    !       j = jloc + xstart(2) - 1
    !       do i = 1, xsize(1)
    !         ub = ub + ux1(i,jloc,k) / ppy(j)
    !       enddo
    !    enddo
    ! enddo

    ! ub = ub * coeff

    ! call MPI_ALLREDUCE(ub,uball,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    
    if(mhd_equation) then
      Je=del_cross_prod(Bm)/Rem
    else
      Je=solve_mhd_potential_poisson(ux1,uy1,uz1)
    endif
    !
    do k = 1, xsize(3)
    do j = 1, xsize(2)
    do i = 1, xsize(1)
      !
      ! Ebar(1)=0.d0
      ! Ebar(2)=0.d0
      ! Ebar(3)=-Bm(i,j,k,2)*uball
      !
      ! elecur(:)=Je(i,j,k,:)+Ebar
      ! elecur(:)=Je(i,j,k,:)
      !
      ! Je(i,j,k,:)
      !
      eforce=cross_product(Je(i,j,k,:),Bm(i,j,k,:))*stuart
      !
      dux1(i,j,k) = dux1(i,j,k)+eforce(1)
      duy1(i,j,k) = duy1(i,j,k)+eforce(2)
      duz1(i,j,k) = duz1(i,j,k)+eforce(3)
      
      ! if(i==1 .and. k==1) then
      !   write(18,*)yy(j),Bm(i,j,k,2),Je(i,j,k,3),eforce(1)
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
  !| this subroutine is used to calculate the result of ∇x.            |
  !+-------------------------------------------------------------------+
  !| change record                                                     |
  !| -------------                                                     |
  !| 28-Oct-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  function del_cross_prod(phi) result(delphi)
    !
    real(mytype) :: delphi(xsize(1),xsize(2),xsize(3),3)
    real(mytype),intent(in) :: phi(xsize(1),xsize(2),xsize(3),3)
    !
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),3) :: dphi
    !
    dphi=grad_vmesh(phi(:,:,:,1))
    !
    delphi(:,:,:,2)= dphi(:,:,:,3)
    delphi(:,:,:,3)=-dphi(:,:,:,2)
    !
    dphi=grad_vmesh(phi(:,:,:,2))
    !
    delphi(:,:,:,1)=-dphi(:,:,:,3)
    delphi(:,:,:,3)= delphi(:,:,:,3) + dphi(:,:,:,1)
    !
    dphi=grad_vmesh(phi(:,:,:,3))
    !
    delphi(:,:,:,1)= delphi(:,:,:,1) + dphi(:,:,:,2)
    delphi(:,:,:,2)= delphi(:,:,:,2) - dphi(:,:,:,1)
    !
    return
    !
  end function del_cross_prod
  !+-------------------------------------------------------------------+
  !| The end of the function del_cross_prod.                           |
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
  function solve_mhd_potential_poisson(ux1,uy1,uz1) result(jcurrent)

    use decomp_2d, only : mytype, xsize, zsize, ph1, nrank
    use decomp_2d_poisson, only : poisson
    use var, only : nzmsize,dv3
    use param, only : ntime, nrhotime, npress,ilmn, ivarcoeff, zero, one 
    use navier,only : gradp

    implicit none

    !! inputs
    real(mytype),dimension(xsize(1), xsize(2), xsize(3)),intent(in) :: ux1, uy1, uz1
    real(mytype),dimension(xsize(1), xsize(2), xsize(3),1:3) :: jcurrent
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
      ucB(i,j,k,:) =cross_product(var1,Bm(i,j,k,:))
      !
    enddo
    enddo
    enddo
    !
    rhs=divergence_scalar(ucB,nlock)
    !
    converged=.false.
    !
    do while(.not.converged)
      !
      call poisson(rhs)
      !
      CALL gradp(jcurrent(:,:,:,1),jcurrent(:,:,:,2),jcurrent(:,:,:,3),rhs)
      !
      converged=.true.
    enddo
    !
    do k = 1, xsize(3)
    do j = 1, xsize(2)
    do i = 1, xsize(1)
      !
      jcurrent(i,j,k,:)=-jcurrent(i,j,k,:)+ucB(i,j,k,:)
      !
    enddo
    enddo
    enddo
    !
    return
    !
    ! div3=divergence_scalar(Je,nlock)
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
  end function solve_mhd_potential_poisson
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solve_mhd_potential_poisson.                      |
  !+-------------------------------------------------------------------+
  !
  subroutine calculate_mhd_transeq_rhs(ux1,uy1,uz1)
    !
    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    !
    call mhd_rhs_eq(dBm(:,:,:,:,1),Bm,ux1,uy1,uz1)
    !
  end subroutine calculate_mhd_transeq_rhs
  !
  subroutine test_magnetic
    !
    use decomp_2d, only : mytype, xsize, zsize, ph1, nrank
    use var,       only : nzmsize,itime,ilist,ifirst,ilast,numscalar, dv3
    use navier, only : divergence
    use param, only : ntime,nrhotime
    !
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: div3


    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar)  :: phi1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    real(mytype), dimension(zsize(1), zsize(2), zsize(3))  :: divu3

    !
    real(mytype) :: maxb,meanb
    integer :: nlock
    !
    nlock=2
    !
    if ((mod(itime,ilist)==0 .or. itime == ifirst .or. itime == ilast)) then
      !
      call divergence(div3,rho1,Bm(:,:,:,1),Bm(:,:,:,2),Bm(:,:,:,3),ep1,drho1,divu3,2,identifier='B')
      !
      !
    endif
    !
  end subroutine test_magnetic
  !
  subroutine mhd_rhs_eq(dB,B,ux1,uy1,uz1)

    use param
    use variables
    use decomp_2d
    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,mu1,mu2,mu3
    use var, only : ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
    use var, only : ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    use var, only : sgsx1,sgsy1,sgsz1
    use var, only : FTx, FTy, FTz, Fdiscx, Fdiscy, Fdiscz
    use ibm_param, only : ubcx,ubcy,ubcz
    use mpi

    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),1:3) :: B,dB
    
    integer :: i,j,k,is
    !

    ! local 
    real(mytype) :: rrem
    real(mytype), save, allocatable, dimension(:,:,:) :: tx1,ty1,tz1,tx2,ty2,tz2, &
                                                         tx3,ty3,tz3,bx2,by2,bz2, &
                                                         bx3,by3,bz3
    logical,save :: firstcal=.true.

    if(firstcal) then
      !
      call alloc_x(tx1)
      tx1 = zero
      call alloc_x(ty1)
      ty1 = zero
      call alloc_x(tz1)
      tz1 = zero
  
      call alloc_y(tx2)
      tx2=zero
      call alloc_y(ty2)
      ty2=zero
      call alloc_y(tz2)
      tz2=zero
  
      call alloc_y(bx2)
      bx2=zero
      call alloc_y(by2)
      by2=zero
      call alloc_y(bz2)
      bz2=zero

      call alloc_z(tx3)
      tx3=zero
      call alloc_z(ty3)
      ty3=zero
      call alloc_z(tz3)
      tz3=zero
      !
      call alloc_z(bx3)
      bx3=zero
      call alloc_z(by3)
      by3=zero
      call alloc_z(bz3)
      bz3=zero
      !
      firstcal=.false.
    endif

    rrem=1.d0/Rem

    !WORK X-PENCILS
    tb1(:,:,:) = ux1(:,:,:) * B(:,:,:,3) - B(:,:,:,1) * uz1(:,:,:)
    tc1(:,:,:) = ux1(:,:,:) * B(:,:,:,2) - B(:,:,:,1) * uy1(:,:,:)

    call derxBy (te1,tb1,di1,sx,ffxB(:,2),fsxB(:,2),fwxB(:,2),xsize(1),xsize(2),xsize(3),0,ubcx*ubcy)
    call derxBz (tf1,tc1,di1,sx,ffxB(:,3),fsxB(:,3),fwxB(:,3),xsize(1),xsize(2),xsize(3),0,ubcx*ubcz)

    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
    
    call transpose_x_to_y(B(:,:,:,1),bx2)
    call transpose_x_to_y(B(:,:,:,2),by2)
    call transpose_x_to_y(B(:,:,:,3),bz2)


    !WORK Y-PENCILS
    ta2(:,:,:) = uy2(:,:,:) * Bz2(:,:,:) - by2 * uz2(:,:,:) 
    tc2(:,:,:) = ux2(:,:,:) * by2(:,:,:) - bx2 * uy2(:,:,:)

    call deryBx (td2,ta2,di2,sy,ffyB(:,1),  fsyB(:,1), fwyB(:,1),ppy,ysize(1),ysize(2),ysize(3),0,ubcx*ubcy)
    call deryBz (tf2,tc2,di2,sy,ffyB(:,3),  fsyB(:,3), fwyB(:,3),ppy,ysize(1),ysize(2),ysize(3),0,ubcz*ubcy)


    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)

    call transpose_y_to_z(bx2,bx3)
    call transpose_y_to_z(by2,by3)
    call transpose_y_to_z(bz2,bz3)

    !WORK Z-PENCILS

    ta3(:,:,:) =  uy3(:,:,:)*bz3(:,:,:) - by3(:,:,:)*uz3(:,:,:)
    tb3(:,:,:) =  ux3(:,:,:)*bz3(:,:,:) - bx3(:,:,:)*uz3(:,:,:)
    tc3(:,:,:) =  ux3(:,:,:)*by3(:,:,:) - bx3(:,:,:)*uy3(:,:,:)


    call derzBx (td3,ta3,di3,sz,ffzB(:,1),fszB(:,1),fwzB(:,1),zsize(1),zsize(2),zsize(3),0,ubcx*ubcz)
    call derzBy (te3,tb3,di3,sz,ffzB(:,2),fszB(:,2),fwzB(:,2),zsize(1),zsize(2),zsize(3),0,ubcy*ubcz)
    call derzBz (tf3,tc3,di3,sz,ffzpB(:,3),fszpB(:,3),fwzpB(:,3),zsize(1),zsize(2),zsize(3),1,ubcz*ubcz)

    !DIFFUSIVE TERMS IN Z
    call derzzBx (ta3,bx3,di3,sz,sfzpB(:,1),sszpB(:,1),swzpB(:,1),zsize(1),zsize(2),zsize(3),1,ubcx)
    call derzzBy (tb3,by3,di3,sz,sfzpB(:,2),sszpB(:,2),swzpB(:,2),zsize(1),zsize(2),zsize(3),1,ubcy)
    call derzzBz (tc3,bz3,di3,sz,sfzB(:,3) ,sszB(:,3) ,swzB(:,3) ,zsize(1),zsize(2),zsize(3),0,ubcz)

    ! Add convective and diffusive terms of z-pencil (half for skew-symmetric)

    td3(:,:,:) = rrem*ta3(:,:,:) - te3(:,:,:)
    te3(:,:,:) = rrem*tb3(:,:,:) + te3(:,:,:)
    tf3(:,:,:) = rrem*tc3(:,:,:)


    !WORK Y-PENCILS
    call transpose_z_to_y(td3,ta2)
    call transpose_z_to_y(te3,tb2)
    call transpose_z_to_y(tf3,tc2)

    ta2(:,:,:) = ta2(:,:,:) + tf2(:,:,:)
    tb2(:,:,:) = tb2(:,:,:) 
    tc2(:,:,:) = tc2(:,:,:) - td2(:,:,:)


    !DIFFUSIVE TERMS IN Y
    if (iimplicit.le.0) then
       !-->for ux
       call deryyBx (td2,bx2,di2,sy,sfypB(:,1),ssypB(:,1),swypB(:,1),ysize(1),ysize(2),ysize(3),1,ubcx)
       if (istret.ne.0) then
          call deryBx (te2,bx2,di2,sy,ffypB(:,1),fsypB(:,1),fwypB(:,1),ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
          do k = 1,ysize(3)
             do j = 1,ysize(2)
                do i = 1,ysize(1)
                   td2(i,j,k) = td2(i,j,k)*pp2y(j)-pp4y(j)*te2(i,j,k)
                enddo
             enddo
          enddo
       endif

       !-->for uy
    call deryyBy (te2,by2,di2,sy,sfyB(:,2),ssyB(:,2),swyB(:,2),ysize(1),ysize(2),ysize(3),0,ubcy)
       if (istret.ne.0) then
          call deryBy (tf2,by2,di2,sy,ffyB(:,2),fsyB(:,2),fwyB(:,2),ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
          do k = 1,ysize(3)
             do j = 1,ysize(2)
                do i = 1,ysize(1)
                   te2(i,j,k) = te2(i,j,k)*pp2y(j)-pp4y(j)*tf2(i,j,k)
                enddo
             enddo
          enddo
       endif

       !-->for uz
    call deryyBz (tf2,bz2,di2,sy,sfypB(:,3),ssypB(:,3),swypB(:,3),ysize(1),ysize(2),ysize(3),1,ubcz)
       if (istret.ne.0) then
          call deryBz (tj2,bz2,di2,sy,ffypB(:,3),fsypB(:,3),fwypB(:,3),ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
          do k = 1,ysize(3)
             do j = 1,ysize(2)
                do i = 1,ysize(1)
                   tf2(i,j,k) = tf2(i,j,k)*pp2y(j)-pp4y(j)*tj2(i,j,k)
                enddo
             enddo
          enddo
       endif
    else ! (semi)implicit Y diffusion
       if (istret.ne.0) then

          !-->for ux
          call deryBx (te2,bx2,di2,sy,ffypB(:,1),fsypB(:,1),fwypB(:,1),ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
          do k=1,ysize(3)
             do j=1,ysize(2)
                do i=1,ysize(1)
                   td2(i,j,k)=-pp4y(j)*te2(i,j,k)
                enddo
             enddo
          enddo
          !-->for uy
          call deryBy (tf2,by2,di2,sy,ffyB(:,2),fsyB(:,2),fwyB(:,2),ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
          do k=1,ysize(3)
             do j=1,ysize(2)
                do i=1,ysize(1)
                   te2(i,j,k)=-pp4y(j)*tf2(i,j,k)
                enddo
             enddo
          enddo
          !-->for uz
          call deryBz (tj2,bz2,di2,sy,ffypB(:,3),fsypB(:,3),fwypB(:,3),ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
          do k=1,ysize(3)
             do j=1,ysize(2)
                do i=1,ysize(1)
                   tf2(i,j,k)=-pp4y(j)*tj2(i,j,k)
                enddo
             enddo
          enddo

       else
       
          td2(:,:,:) = zero
          te2(:,:,:) = zero
          tf2(:,:,:) = zero
          
       endif
    endif

    ! Add diffusive terms of y-pencil to convective and diffusive terms of y- and z-pencil
    ta2(:,:,:) = rrem*td2(:,:,:) + ta2(:,:,:)
    tb2(:,:,:) = rrem*te2(:,:,:) + tb2(:,:,:)
    tc2(:,:,:) = rrem*tf2(:,:,:) + tc2(:,:,:)

    !WORK X-PENCILS
    call transpose_y_to_x(ta2,ta1)
    call transpose_y_to_x(tb2,tb1)
    call transpose_y_to_x(tc2,tc1) !diff+conv. terms

    !DIFFUSIVE TERMS IN X
    call derxxBx (td1,B(:,:,:,1),di1,sx,sfxB(:,1) ,ssxB(:,1) ,swxB(:,1) ,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derxxBy (te1,B(:,:,:,2),di1,sx,sfxpB(:,2),ssxpB(:,2),swxpB(:,2),xsize(1),xsize(2),xsize(3),1,ubcy)
    call derxxBz (tf1,B(:,:,:,3),di1,sx,sfxpB(:,3),ssxpB(:,3),swxpB(:,3),xsize(1),xsize(2),xsize(3),1,ubcz)

    !FINAL SUM: DIFF TERMS + CONV TERMS
    dB(:,:,:,1) = rrem * td1(:,:,:) + ta1(:,:,:) 
    dB(:,:,:,2) = rrem * te1(:,:,:) + tb1(:,:,:) -tf1(:,:,:)
    dB(:,:,:,3) = rrem * te1(:,:,:) + tc1(:,:,:) +te1(:,:,:)

    return

  end subroutine mhd_rhs_eq
  !
  subroutine solve_poisson_mhd
    !
    use decomp_2d, only : mytype, xsize, zsize, ph1, nrank
    use decomp_2d_poisson, only : poisson
    use var, only : nzmsize,dv3
    use param, only : ntime, nrhotime, npress,ilmn, ivarcoeff, zero, one 
    use navier,only : gradp

    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: phib

    real(mytype),dimension(xsize(1),xsize(2),xsize(3),1:3) :: dphib

    integer :: i,j,k,nlock,poissiter
    !
    nlock=1 !! Corresponds to computing div(u*)
    !
    do poissiter = 1, 1
      phib=divergence_scalar(Bm,nlock) !todo: this will have incorrect BCs?
      call poisson(phib)
      CALL gradp(dphib(:,:,:,1),dphib(:,:,:,2),dphib(:,:,:,3),phib)
      Bm=Bm-dphib
    enddo
    !
    sync_Bm_needed=.true.
    !
  end subroutine solve_poisson_mhd
  !
  ! subroutine solve_poisson_mhd2(rho1, ep1, drho1, divu3)

  !   USE decomp_2d, ONLY : mytype, xsize, zsize, ph1, nrank, real_type
  !   USE decomp_2d_poisson, ONLY : poisson
  !   USE var, ONLY : nzmsize,dv3
  !   USE param, ONLY : ntime, nrhotime, npress,ilmn, ivarcoeff, zero, one 
  !   USE mpi
  !   use navier,only : gradp,velocity_to_momentum,divergence
    
  !   implicit none

  !   !! Inputs
  !   REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ep1
  !   REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime), INTENT(IN) :: rho1
  !   REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime), INTENT(IN) :: drho1
  !   REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)), INTENT(IN) :: divu3

  !   !! Outputs
  !   REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
  !   REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1

  !   !! Locals
  !   INTEGER :: nlock, poissiter
  !   LOGICAL :: converged
  !   REAL(mytype) :: atol, rtol, rho0, divup3norm


  !   nlock = 1 !! Corresponds to computing div(u*)
  !   converged = .FALSE.
  !   poissiter = 0
  !   rho0 = one



  !   CALL divergence(pp3(:,:,:,1),rho1,Bm(:,:,:,1),Bm(:,:,:,2),Bm(:,:,:,3),ep1,drho1,divu3,nlock,identifier='B')
  !   !
  !   IF (ilmn.AND.ivarcoeff) THEN
  !      dv3(:,:,:) = pp3(:,:,:,1)
  !   ENDIF

  !   do while(.not.converged)

  !      IF (.NOT.converged) THEN
  !         CALL poisson(pp3(:,:,:,1))


  !         !! Need to update pressure gradient here for varcoeff
  !         CALL gradp(px1,py1,pz1,pp3(:,:,:,1))

         

  !         IF ((.NOT.ilmn).OR.(.NOT.ivarcoeff)) THEN
  !            !! Once-through solver
  !            !! - Incompressible flow
  !            !! - LMN - constant-coefficient solver
  !            converged = .TRUE.
  !         ENDIF
  !      ENDIF

  !      poissiter = poissiter + 1
  !   enddo


  !   Bm(:,:,:,1)=Bm(:,:,:,1)-px1
  !   Bm(:,:,:,2)=Bm(:,:,:,2)-py1
  !   Bm(:,:,:,3)=Bm(:,:,:,3)-pz1

  !   sync_Bm_needed=.true.
  !   !
  ! end subroutine solve_poisson_mhd2
  !
  subroutine  div_check(divec,divmax,divmean)
    !
    USE decomp_2d
    !
    real(mytype),intent(in) :: divec(:,:,:)
    real(mytype),intent(out) :: divmax,divmean
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
    divmax=tmax
    divmean=tmoy

    ! if (nrank == 0) then
    !    write(*,*) 'DIV max mean=',real(tmax,mytype),real(tmoy/real(nproc),mytype)
    ! endif
    !
  end subroutine
  !
  subroutine  divergence_vmesh_check(vec,divec)
    !
    USE MPI
    USE decomp_2d
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
       write(*,*) 'DIV max mean=',real(tmax,mytype),real(tmoy/real(nproc),mytype)
    endif
    !
  end subroutine  divergence_vmesh_check
  !+-------------------------------------------------------------------+
  !| The end of the subroutine divergence_vmesh_check.                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is used to check the velocity at wall             |
  !+-------------------------------------------------------------------+
  !| change record                                                     |
  !| -------------                                                     |
  !| 27-Jan-2023  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine boundary_velocity_check(ux,uy,uz,note)
    !
    USE decomp_2d
    use var,     only : nx,ny,nz
    !
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) ::   &
                                                             ux,uy,uz
    character(len=*),intent(in),optional :: note
    !
    real(8) :: ux_m,uy_m,uz_m,ux_f,uy_f,uz_f
    integer :: i,j,k
    !
    ux_m=0.d0; uy_m=0.d0; uz_m=0.d0
    ux_f=0.d0; uy_f=0.d0; uz_f=0.d0
    !
    if (xstart(2)==1) then
      do k=1,xsize(3)
      do i=1,xsize(1)
        ux_m=ux_m+ux(i,1,k)
        ux_f=ux_f+ux(i,1,k)**2
        !
        uy_m=uy_m+uy(i,1,k)
        uy_f=uy_f+uy(i,1,k)**2
        !
        uz_m=uz_m+uz(i,1,k)
        uz_f=uz_f+uz(i,1,k)**2
      enddo
      enddo
    endif
    !
    if (xend(2)==ny) then
      do k=1,xsize(3)
      do i=1,xsize(1)
        ux_m=ux_m+ux(i,xsize(2),k)
        ux_f=ux_f+ux(i,xsize(2),k)**2
        !
        uy_m=uy_m+uy(i,xsize(2),k)
        uy_f=uy_f+uy(i,xsize(2),k)**2
        !
        uz_m=uz_m+uz(i,xsize(2),k)
        uz_f=uz_f+uz(i,xsize(2),k)**2
      enddo
      enddo
    endif
    !
    ux_m=psum(ux_m)/(nx*nz)
    ux_f=psum(ux_f)/(nx*nz)
    uy_m=psum(uy_m)/(nx*nz)
    uy_f=psum(uy_f)/(nx*nz)
    uz_m=psum(uz_m)/(nx*nz)
    uz_f=psum(uz_f)/(nx*nz)
    !
    if (nrank == 0) then
      !
      if(present(note)) then
        write(*,*)' ** note ',note
      endif
      !
      write(*,*)' ** u_m ', ux_m,uy_m,uz_m
      write(*,*)' ** u_f ', ux_f,uy_f,uz_f
    endif
    !
  end subroutine boundary_velocity_check
  !+-------------------------------------------------------------------+
  !| The end of the subroutine boundary_velocity_check.                |
  !+-------------------------------------------------------------------+
  !

  !!############################################################################
  !subroutine DIVERGENCe
  !Calculation of div u* for nlock=1 and of div u^{n+1} for nlock=2
  ! input :  vec (on velocity mesh)
  ! output : pp3 (on pressure mesh)
  !written by SL 2018
  !############################################################################
  function divergence_scalar(vec,nlock) result(pp3)

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

    tmax=-1609._mytype
    tmoy=zero
    do k=1,nzmsize
       do j=ph1%zst(2),ph1%zen(2)
          do i=ph1%zst(1),ph1%zen(1)
             if (pp3(i,j,k).gt.tmax) tmax=pp3(i,j,k)
             tmoy=tmoy+abs(pp3(i,j,k))
          enddo
       enddo
    enddo
    tmoy=tmoy/nvect3

    call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    if ((nrank == 0) .and. (nlock > 0).and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
       if (nlock == 2) then
          write(*,*) 'DIV B  max mean=',real(tmax1,mytype),real(tmoy1/real(nproc),mytype)
       else
          write(*,*) 'DIV B* max mean=',real(tmax1,mytype),real(tmoy1/real(nproc),mytype)
       endif
    endif

    return
    !
  end function divergence_scalar
  !
end module mhd
!+---------------------------------------------------------------------+
! the end of the module mhd                                            |
!+---------------------------------------------------------------------+
