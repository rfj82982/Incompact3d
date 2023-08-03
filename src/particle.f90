!+---------------------------------------------------------------------+
!| This module contains variables, arraies and subroutines related to  |
!| particle tracking.                                                  |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 15-11-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module partack
  !
  use decomp_2d, only : mytype,nrank,nproc
  use mptool
  use constants, only : pi
  !
  implicit none
  !
  interface mclean
    module procedure mclean_mytype
    module procedure mclean_particle
  end interface mclean 
  !
  interface msize
    module procedure size_particle
    module procedure size_integer
  end interface msize
  !
  interface pa2a
    module procedure pa2a_particle
  end interface pa2a
  !
  interface mextend
     module procedure extend_particle
  end interface mextend
  !
  ! particles
  type partype
    !
    real(mytype) :: rho,mas,dim,re,vdiff,qe
    real(mytype) :: x(3),v(3),vf(3),f(3),b(3)
    real(mytype),allocatable,dimension(:,:) :: dx,dv
    integer :: id,rankinn,rank2go
    logical :: swap,new
    !+------------------+------------------------------------------+
    !|              rho | density                                  |
    !|              mas | mass                                     |
    !|              dim | dimeter                                  |
    !|               Re | particle Reynolds number.                |
    !|            vdiff | velocity difference between fluid and    |
    !|                  | particle                                 |
    !|               qe | Electric charge                          |
    !|                x | spatial coordinates of particle          |
    !|                v | velocity of particle                     |
    !|               vf | velocity of fluids                       |
    !|                f | force                                    |
    !|                b | magnetic field                           |
    !|               dx | gradient of x,y,z to time, used for      |
    !|                  | temporal integration.                    |
    !|               dv | gradient of u,v,w to time, used for      |
    !|                  | temporal integration.                    |
    !|               id | the identification of particle           |
    !|          rankinn | the mpi rank which the particle is in    |
    !|          rank2go | the mpi rank which the particle will be  |
    !+------------------+------------------------------------------+
    !
    contains
    !
    procedure :: init  => init_one_particle
    procedure :: reset => reset_one_particle
    procedure :: rep   => particle_reynolds_cal
    procedure :: force => particle_force_cal
    !
  end type partype
  !
  type(partype),allocatable,target :: particle(:)
  !
  logical :: lpartack,lorentzforce,lfluidforce
  integer :: numparticle,ipartiout,ipartiadd,particle_file_numb
  real(mytype) :: partirange(6)
  integer :: numpartix(3)
  real(mytype),allocatable,dimension(:) :: lxmin,lxmax,lymin,lymax,lzmin,lzmax
  real(mytype),allocatable,dimension(:) :: xpa,ypa,zpa
  real(mytype),allocatable,dimension(:) :: ux_pa,uy_pa,uz_pa
  real(mytype) :: sub_time_step,particletime
  real(mytype) :: part_time,part_comm_time,part_vel_time,part_dmck_time,a2a_time, &
                  count_time,data_pack_time,data_unpack_time,mpi_comm_time,       &
                  table_share_time,h5io_time
  !+------------------+--------------------------------------------+
  !|         lpartack | switch of particel tracking                |
  !|      numparticle | number of particles in the domain          |
  !|        ipartiout | frequency of output particles              |
  !|        ipartiadd | frequency of add new particles             |
  !|       partirange | the domain where the particles are injucted|
  !|        numpartix | the matrix of particle number              |
  !|            ux_pa |                                            |
  !|            uy_pa |                                            |
  !|            uz_pa | velocity of particles                      |
  !+------------------+--------------------------------------------+
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to init a particle.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine init_one_particle(pa)
    !
    use param, only : ntime
    !
    class(partype),target :: pa
    !
    pa%swap=.false.
    pa%new =.true.
    !
    pa%rankinn=nrank
    pa%rank2go=nrank
    !
    pa%x=0.0; pa%v=0.0
    !
    allocate(pa%dx(1:3,2:ntime),pa%dv(1:3,ntime))
    !
    pa%dx=0.0
    pa%dv=0.0
    !
    pa%dim=1.d-1
    pa%rho=10._mytype
    !
    pa%qe=1.d-5
    !
  end subroutine init_one_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine initmesg.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to reset a particle.                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-Jun-2022  | Created by J. Fang @ Imperial College              |
  !+-------------------------------------------------------------------+
  subroutine reset_one_particle(pa)
    !
    use param, only : ntime
    !
    class(partype),target :: pa
    !
    pa%swap=.false.
    pa%new =.true.
    !
    pa%rankinn=nrank
    pa%rank2go=nrank
    !
    pa%x=0.0; pa%v=0.0
    !
    pa%dx=0.0
    pa%dv=0.0
    !
  end subroutine reset_one_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine reset_one_particle.                     |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate the Reynolds number of       |
  !| particles.                                                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Nov-2022  | Created by J. Fang @ Imperial College              |
  !+-------------------------------------------------------------------+
  subroutine particle_reynolds_cal(pa)
    !
    use param, only : re
    !
    class(partype),target :: pa
    !
    pa%vdiff = sqrt( (pa%vf(1)-pa%v(1))**2 + &
                     (pa%vf(2)-pa%v(2))**2 + &
                     (pa%vf(3)-pa%v(3))**2 )
    !
    pa%re = pa%dim*pa%vdiff*re
    !
  end subroutine particle_reynolds_cal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_reynolds_cal.                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate the foce acting on a particle|                                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Nov-2022  | Created by J. Fang @ Imperial College              |
  !+-------------------------------------------------------------------+
  subroutine particle_force_cal(pa)
    !
    use param, only : re
    !
    class(partype),target :: pa
    !
    real(mytype) :: varc
    ! 
    varc=18.d0/(pa%rho*pa%dim**2*re)*(1.d0+0.15d0*pa%re**0.687d0)
    !
    pa%f(:) = varc*(pa%vf(:)-pa%v(:))
    !
    print*,pa%f(:),varc,pa%vf(:)-pa%v(:)
    ! !
    ! call mpistop
    !
  end subroutine particle_force_cal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_reynolds_cal.                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to report time cost for particles.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine partcle_report
    !
    integer :: ttal_particle
    !
    ttal_particle=psum(numparticle)
    !
    if(nrank==0) then
      write(*,*) 'Total number of particles:',ttal_particle
      write(*,*) 'Total time for particles :',real(part_time,4)
      write(*,*) '      time particles vel :',real(part_vel_time,4)
      write(*,*) '      time domain search :',real(part_dmck_time,4)
      write(*,*) '      time partical_swap :',real(part_comm_time,4)
      write(*,*) '           alltoall comm :',real(a2a_time,4)
      write(*,*) '           counting time :',real(count_time,4)
      write(*,*) '           table shareing:',real(table_share_time,4)
      write(*,*) '           data packing  :',real(data_pack_time,4)
      write(*,*) '           MPI Alltoall  :',real(mpi_comm_time,4)
      write(*,*) '           data unpacking:',real(data_unpack_time,4)
      write(*,*) '                  hdf5 io:',real(h5io_time,4)
    endif
    !
  end subroutine partcle_report
  !+-------------------------------------------------------------------+
  !| The end of the subroutine partcle_report.                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate particle array.                    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_gen(particle_new,particle_size)
    !
    ! arguments
    type(partype),intent(out),allocatable :: particle_new(:)
    integer,intent(out) :: particle_size
    !
    ! local data
    integer :: p,i,j,k,max_part_size
    real(mytype) :: dx,dy,dz,x,y,z
    !
    p=0
    !
    max_part_size=numpartix(1)*numpartix(2)*numpartix(3)
    allocate(particle_new(1:max_part_size))
    !
    do k=1,numpartix(3)
    do j=1,numpartix(2)
    do i=1,numpartix(1)
      !
      dx=(partirange(2)-partirange(1))/real(numpartix(1),mytype)
      dy=(partirange(4)-partirange(3))/real(numpartix(2),mytype)
      dz=(partirange(6)-partirange(5))/real(numpartix(3),mytype)
      !
      x=dx*real(i,mytype)+partirange(1)
      y=dy*real(j,mytype)+partirange(3)
      z=dz*real(k,mytype)+partirange(5)
      !
      ! print*,x,y,z
      ! print*,(y>=lymin(nrank) .and. y<=lymax(nrank))
      !
      if( x>=lxmin(nrank) .and. x<=lxmax(nrank) .and. &
          y>=lymin(nrank) .and. y<=lymax(nrank) .and. &
          z>=lzmin(nrank) .and. z<=lzmax(nrank) ) then
        !
        p=p+1
        !
        call particle_new(p)%init()
        !
        particle_new(p)%x(1)=x
        particle_new(p)%x(2)=y
        particle_new(p)%x(3)=z
        !
        particle_new(p)%v(1)   =y !0.d0
        particle_new(p)%v(2)   =0.d0
        particle_new(p)%v(3)   =0.d0
        !
        particle_new(p)%new=.false.
        !
      endif
      !
    enddo
    enddo
    enddo
    !
    call mclean(particle_new,p)
    !
    particle_size=p
    !
  end subroutine particle_gen
  !+-------------------------------------------------------------------+
  !| The end of the subroutine particle_gen.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to initilise particle positions.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine init_particle
    !
    use param,     only : xlx,yly,zlz,irestart,t,dt
    use var,       only : itime
    !
    ! local data
    integer :: i,j,k,p
    real(mytype) :: dx,dy,dz
    !
    if(irestart==0) then
      !
      call particle_gen(particle,numparticle)
      !
      particle_file_numb=0
      !
      call h5write_particle()
      !
    else
      ! call h5read_particle(particle,numparticle)
      call particle_gen(particle,numparticle)
      !
      particle_file_numb=0
      !
      call h5write_particle()
      !
    endif
    !
    ! call partical_domain_check('bc_tgv')
    ! call partical_domain_check('bc_channel')
    call partical_domain_check('out_disappear')
    !
    !
    ! call partical_swap
    !
    ! call write_particle()
    !
    part_time=0.d0
    part_comm_time=0.d0
    part_vel_time=0.d0
    part_dmck_time=0.d0
    a2a_time=0.d0
    count_time=0.d0
    table_share_time=0.d0
    data_pack_time=0.d0
    data_unpack_time=0.d0
    mpi_comm_time=0.d0
    h5io_time=0.d0
    !
    particletime=t
    sub_time_step=dt
    ! sub_time_step=0.1d0*dt
    !
    ! lorentzforce=.true.
    !
  end subroutine init_particle
  !+-------------------------------------------------------------------+
  ! The end of the subroutine init_particle                            |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to add more particles to the domain.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine partile_inject
    !
    use param,     only : xlx,yly,zlz
    !
    !
    ! local data
    type(partype),allocatable :: particle_new(:)
    integer :: num_new_particle,n
    !
    call particle_gen(particle_new,num_new_particle)
    !
    call particle_add(particle,particle_new,n)
    !
    numparticle=numparticle+n
    !
    ! call partical_domain_check('out_disappear')
    call partical_domain_check('out_disappear')
    !
    call partical_swap
    !
  end subroutine partile_inject
  !+-------------------------------------------------------------------+
  ! The end of the subroutine partile_add                              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to calcualte the size and range of local domain|
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine local_domain_size
    !
    use param,     only : dx,dy,dz,istret
    use variables, only : yp,ny
    use decomp_2d, only : xstart,xend
    use actuator_line_model_utils
    !
    integer :: nyr,nzr,jrank
    !
    allocate( lxmin(0:nproc-1),lxmax(0:nproc-1),  &
              lymin(0:nproc-1),lymax(0:nproc-1),  &
              lzmin(0:nproc-1),lzmax(0:nproc-1)   )
    !
    lxmin=0.d0; lxmax=0.d0
    lymin=0.d0; lymax=0.d0
    lzmin=0.d0; lzmax=0.d0
    !

    lxmin(nrank)=real(0,mytype)*dx
    lxmax(nrank)=real(xend(1),mytype)*dx
    !
    if (istret==0) then
      lymin(nrank)=real(xstart(2)-1,mytype)*dy
      lymax(nrank)=real(xend(2),mytype)*dy
    else
      lymin(nrank)=yp(xstart(2))
      nyr=min(ny,xend(2)+1)
      lymax(nrank)=yp(nyr)
    endif
    !
    lzmin(nrank)=real((xstart(3)-1),mytype)*dz
    nzr=xend(3)
    lzmax(nrank)=real(nzr,mytype)*dz
    !
    lxmin=psum(lxmin); lxmax=psum(lxmax)
    lymin=psum(lymin); lymax=psum(lymax)
    lzmin=psum(lzmin); lzmax=psum(lzmax)
    !
    ! if(nrank==5) then
    !   do jrank=0,nproc-1
    !     print*,nrank,jrank,lxmin(jrank),lxmax(jrank),lymin(jrank),lymax(jrank),lzmin(jrank),lzmax(jrank)
    !   enddo
    ! endif
    !
    ! call mpistop
    !
  end subroutine local_domain_size
  !+-------------------------------------------------------------------+
  ! The end of the subroutine local_domain_size                        |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate the force on particles.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-11-2021  | Created by J. Fang                                  |
  !| 01-04-2023  | Add lorentz force by J. Fang                        |
  !+-------------------------------------------------------------------+
  subroutine particle_force(ux1,uy1,uz1)
    !
    use decomp_2d, only : xsize
    use var,       only : itime
    use param,     only : re,t,itr
    use mhd,       only : Bm,stuart
    !
    ! arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    !
    ! local data
    integer :: psize,jpart,npart
    type(partype),pointer :: pa
    real(mytype) :: varc,vard,varq(3),maxforce,maxvdiff,vdiff,re_p,bmt(3)
    !
    logical,save :: firstcal=.true.
    !
    if(firstcal) then
      !
      ! if(nrank==0) open(52,file='particle_log')
      !
      if(.not. allocated(bm)) allocate(bm(xsize(1),xsize(2),xsize(3),1:3))
      !
      bm(:,:,:,1)=0.d0
      bm(:,:,:,2)=1.d0
      bm(:,:,:,3)=0.d0
      !
      firstcal=.false.
      !
    endif
    !
    bmt(1)=0.d0
    bmt(2)=1.d0
    bmt(3)=0.d0
    !
    if(lorentzforce) then
      call field_var(ux1,uy1,uz1,Bm)
    else
      call field_var(ux1,uy1,uz1)
    endif
    !
    psize=msize(particle)
    !
    ! print*,nrank,'-',numparticle,'-',psize
    !
    npart=0
    !
    maxforce=0.d0
    maxvdiff=0.d0
    !
    ! print*,nrank,'|',numparticle
    !
    do jpart=1,psize
      !
      pa=>particle(jpart)
      !
      if(pa%new) cycle
      !
      ! ! tracking particles
      ! pa%v=pa%vf
      ! pa%f(:)=0.d0
      !
      pa%f(:) = 0.d0
      !
      if(lfluidforce) then

        vdiff = norm2(pa%vf-pa%v)
        re_p = pa%dim*vdiff*re
        varc=18.d0/(pa%rho*pa%dim**2*re)*(1.d0+0.15d0*re_p**0.687d0)
        pa%f(:) = pa%f(:) + varc*(pa%vf(:)-pa%v(:))

      endif
      !
      if(lorentzforce) then

        vard = 6.d0/(pi*(pa%dim**3)*pa%rho)*stuart*pa%qe

        varq(:) = cross_product(pa%v(:),bmt(:))

        pa%f(:) = pa%f(:) + vard*varq(:)

        if(isnan(pa%f(1)).or. isnan(pa%f(2)).or.isnan(pa%f(3))) then
          print*,' !! error in calculating force !!'
          print*,'  varc:',varc,'  vard:',vard
          print*,'    vf:',pa%vf
          print*,'     v:',pa%v
          print*,'     b:',pa%b
          print*,'    vq:',varq
          stop
        endif

        ! print*,jpart,'|',varc*varq(:)

      endif
      !
      if(lfluidforce .or. lorentzforce) then
        continue
      else
        ! no force acting on particles, purely fluid tracker
        pa%v=pa%vf
      endif
      !
      maxforce=max(maxforce,abs(pa%f(1)),abs(pa%f(2)),abs(pa%f(3)))
      maxvdiff=max(maxvdiff,vdiff)
      ! !
      npart=npart+1
      !
      if(npart==numparticle) exit
      !
    enddo
    !
    !
    if(mod(itime, ipartiout).eq.0) then
      maxforce=pmax(maxforce)
      maxvdiff=pmax(maxvdiff)
      !
      if(nrank==0) then
    !     ! write(52,*)particletime,maxforce,maxvdiff
        print*,' ** particle maxforce:',maxforce,'max vdiff:',maxvdiff
      endif
      !
    endif
    !
  end subroutine particle_force
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_force                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to get the fluids velocity on particles.       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine field_var(ux1,uy1,uz1,B)
    !
    use MPI
    use param,     only : dx,dy,dz,istret,nclx,ncly,nclz,xlx,yly,zlz
    use variables, only : yp,ny,nz
    use decomp_2d, only : xsize,xstart,xend,update_halo
    use actuator_line_model_utils
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),1:3),intent(in),optional :: B
    !
    ! local data
    integer :: jpart,npart,i,j,k,psize
    type(partype),pointer :: pa
    real(mytype) :: x1,y1,z1,x2,y2,z2
    real(mytype) :: test(xsize(1),xsize(2),xsize(3))
    real(mytype),allocatable,dimension(:,:,:) :: ux1_halo,uy1_halo,    &
                                                 uz1_halo,ux1_hal2,    &
                                                 uy1_hal2,uz1_hal2
    real(mytype),allocatable :: b_halo(:,:,:,:)
    !
    real(mytype),allocatable,save :: xx(:),yy(:),zz(:)
    logical,save :: firstcal=.true.
    !
    real(mytype) :: timebeg
    !
    timebeg=ptime()
    !
    if(firstcal) then
      !
      allocate(xx(0:xsize(1)+1),yy(0:xsize(2)+1),zz(0:xsize(3)+1))
      !
      do i=0,xsize(1)+1
        xx(i)=real(i-1,mytype)*dx
      enddo
      do j=0,xsize(2)+1
        !
        if(j+xstart(2)-1>ny) then
          yy(j)=2.d0*yp(ny)-yp(ny-1)
        elseif(j+xstart(2)-1<1) then
          yy(j)=2.d0*yp(1)-yp(2)
        else
          yy(j)=yp(j+xstart(2)-1)
        endif
        !
      enddo
      !
      if(xsize(3)==1) then
        zz(:)=0.d0
      else
        do k=0,xsize(3)+1
          zz(k)=real((k+xstart(3)-2),mytype)*dz
        enddo
      endif
      !
      firstcal=.false.
      !
    endif
    !
    allocate( ux1_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1),        &
              uy1_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1),        &
              uz1_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1) )
    ! 
    call pswap_yz(ux1,ux1_halo)
    call pswap_yz(uy1,uy1_halo)
    call pswap_yz(uz1,uz1_halo)
    !
    if(present(B)) then
      !
      allocate( b_halo(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1,1:3) )
      !
      call pswap_yz(b(:,:,:,1),b_halo(:,:,:,1))
      call pswap_yz(b(:,:,:,2),b_halo(:,:,:,2))
      call pswap_yz(b(:,:,:,3),b_halo(:,:,:,3))
      !
    endif
    !
    psize=msize(particle)
    !
    npart=0
    !
    do jpart=1,psize
      !
      pa=>particle(jpart)
      !
      if(pa%new) cycle
      !
      loopk: do k=1,xsize(3)
        !
        z1=zz(k)
        z2=zz(k+1)
        !
        do j=1,xsize(2)
          y1=yy(j)
          y2=yy(j+1)
          do i=1,xsize(1)
            x1=xx(i)
            x2=xx(i+1)
            !
            if( pa%x(1)>=x1 .and. pa%x(1)<=x2 .and. &
                pa%x(2)>=y1 .and. pa%x(2)<=y2 .and. &
                pa%x(3)>=z1 .and. pa%x(3)<=z2 ) then
              !
              ! locate the particle, do the interpolation
              ! print*,x1,x2,y1,y2,z1,z2
              pa%vf(1)=trilinear_interpolation( x1,y1,z1,            &
                                                x2,y2,z2,            &
                                            pa%x(1),pa%x(2),pa%x(3), &
                                            ux1_halo(i,j,k),     &
                                            ux1_halo(i+1,j,k),   &
                                            ux1_halo(i,j,k+1),   &
                                            ux1_halo(i+1,j,k+1), &
                                            ux1_halo(i,j+1,k),   &
                                            ux1_halo(i+1,j+1,k), &
                                            ux1_halo(i,j+1,k+1), &
                                            ux1_halo(i+1,j+1,k+1))
              !
              if(isnan(pa%vf(1))) then
                print*,' !! error in calculating vf x'
                print*,'        x1:',x1,y1,z1
                print*,'        x2:',x2,y2,z2
                print*,'      pa%x:',pa%x
                print*,'      rank:',nrank,'i,j,k',i,j,k
                print*,'  ux1_halo:',ux1_halo(i:i+1,j:j+1,k:k+1)
                stop
              endif
              !
              pa%vf(2)=trilinear_interpolation( x1,y1,z1,           &
                                               x2,y2,z2,            &
                                            pa%x(1),pa%x(2),pa%x(3),&
                                            uy1_halo(i,j,k),     &
                                            uy1_halo(i+1,j,k),   &
                                            uy1_halo(i,j,k+1),   &
                                            uy1_halo(i+1,j,k+1), &
                                            uy1_halo(i,j+1,k),   &
                                            uy1_halo(i+1,j+1,k), &
                                            uy1_halo(i,j+1,k+1), &
                                            uy1_halo(i+1,j+1,k+1)) 
              !
              if(isnan(pa%vf(2))) then
                print*,' !! error in calculating vf y'
                print*,'        x1:',x1,y1,z1
                print*,'        x2:',x2,y2,z2
                print*,'      pa%x:',pa%x
                print*,'      rank:',nrank,'i,j,k',i,j,k
                print*,'  uy1_halo:',uy1_halo(i:i+1,j:j+1,k:k+1)
                stop
              endif
              !
              pa%vf(3)=trilinear_interpolation( x1,y1,z1,            &
                                                x2,y2,z2,            &
                                            pa%x(1),pa%x(2),pa%x(3), &
                                            uz1_halo(i,j,k),     &
                                            uz1_halo(i+1,j,k),   &
                                            uz1_halo(i,j,k+1),   &
                                            uz1_halo(i+1,j,k+1), &
                                            uz1_halo(i,j+1,k),   &
                                            uz1_halo(i+1,j+1,k), &
                                            uz1_halo(i,j+1,k+1), &
                                            uz1_halo(i+1,j+1,k+1)) 
              !
              if(isnan(pa%vf(3))) then
                print*,' !! error in calculating vf z'
                print*,'        x1:',x1,y1,z1
                print*,'        x2:',x2,y2,z2
                print*,'      pa%x:',pa%x
                print*,'      rank:',nrank,'i,j,k',i,j,k
                print*,'  uz1_halo:',uz1_halo(i:i+1,j:j+1,k:k+1)
                stop
              endif
              !
              !
              if(present(B)) then
                !
                pa%b(1)=trilinear_interpolation( x1,y1,z1,            &
                                                  x2,y2,z2,            &
                                              pa%x(1),pa%x(2),pa%x(3), &
                                              b_halo(i,j,k,1),     &
                                              b_halo(i+1,j,k,1),   &
                                              b_halo(i,j,k+1,1),   &
                                              b_halo(i+1,j,k+1,1), &
                                              b_halo(i,j+1,k,1),   &
                                              b_halo(i+1,j+1,k,1), &
                                              b_halo(i,j+1,k+1,1), &
                                              b_halo(i+1,j+1,k+1,1))
                pa%b(2)=trilinear_interpolation( x1,y1,z1,            &
                                                 x2,y2,z2,            &
                                              pa%x(1),pa%x(2),pa%x(3),&
                                              b_halo(i,j,k,2),        &
                                              b_halo(i+1,j,k,2),   &
                                              b_halo(i,j,k+1,2),   &
                                              b_halo(i+1,j,k+1,2), &
                                              b_halo(i,j+1,k,2),   &
                                              b_halo(i+1,j+1,k,2), &
                                              b_halo(i,j+1,k+1,2), &
                                              b_halo(i+1,j+1,k+1,2) )
                pa%b(3)=trilinear_interpolation( x1,y1,z1,            &
                                                  x2,y2,z2,            &
                                              pa%x(1),pa%x(2),pa%x(3), &
                                              b_halo(i,j,k,3),     &
                                              b_halo(i+1,j,k,3),   &
                                              b_halo(i,j,k+1,3),   &
                                              b_halo(i+1,j,k+1,3), &
                                              b_halo(i,j+1,k,3),   &
                                              b_halo(i+1,j+1,k,3), &
                                              b_halo(i,j+1,k+1,3), &
                                              b_halo(i+1,j+1,k+1,3)) 
              if(isnan(pa%b(1)) .or. isnan(pa%b(2)) .or. isnan(pa%b(3))) then
                print*,' !! error in calculating b'
                print*,'        x1:',x1,y1,z1
                print*,'        x2:',x2,y2,z2
                print*,'      pa%x:',pa%x
                print*,'      rank:',nrank,'i,j,k',i,j,k
                print*,'  b_halo_x:',b_halo(i:i+1,j:j+1,k:k+1,1)
                print*,'  b_halo_y:',b_halo(i:i+1,j:j+1,k:k+1,2)
                print*,'  b_halo_z:',b_halo(i:i+1,j:j+1,k:k+1,3)
                stop
              endif
              !
              !
              endif
              !
              exit loopk
              !
            endif
            !
          enddo
        enddo
      enddo loopk
      !
      npart=npart+1
      !
      if(npart==numparticle) exit
      !
    enddo
    !
    part_vel_time=part_vel_time+ptime()-timebeg
    !
    part_time=part_time+ptime()-timebeg
    !
  end subroutine field_var
  !+-------------------------------------------------------------------+
  ! The end of the subroutine field_var                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to integrate particle coordinates in time.     |
  !+-------------------------------------------------------------------+
  !| only Euler scheme is used for now.                                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine intt_particel(ux1,uy1,uz1,time1)
    !
    use variables 
    use param
    use decomp_2d, only : xsize
    !
    ! arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    real(mytype),intent(in) :: time1
    !
    ! local data 
    integer :: p,psize,jpart,npart,total_num_part
    integer,save :: old_num_part=0
    type(partype),pointer :: pa
    real(mytype) :: timebeg
    real(mytype),save :: time0,tsro
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxp,uyp,uzp
    real(mytype),allocatable,dimension(:,:,:),save :: ux0,uy0,uz0
    real(mytype),allocatable :: xcor(:,:),dxco(:,:,:),vcor(:,:),dvco(:,:,:)
    !
    logical,save :: firstcal=.true.
    !
    if(firstcal) then
      !
      allocate( ux0(xsize(1),xsize(2),xsize(3)), &
                uy0(xsize(1),xsize(2),xsize(3)), &
                uz0(xsize(1),xsize(2),xsize(3)) )
      !
      ux0=ux1
      uy0=uy1
      uz0=uz1
      !
      time0=time1
      !
      tsro=sub_time_step/dt
      !
      particletime=time1
      !
      firstcal=.false.
      !
      return
      !
    endif
    !
    timebeg=ptime()
    !
    do while(particletime<time1)
      !
      particletime=particletime+sub_time_step
      !
      ! print*,'time0=',time0,'time1=',time1,'particletime=',particletime
      !
      uxp=linintp(time0,time1,ux0,ux1,particletime)
      uyp=linintp(time0,time1,uy0,uy1,particletime)
      uzp=linintp(time0,time1,uz0,uz1,particletime)
      !
      call particle_force(uxp,uyp,uzp)
      !
      allocate(xcor(3,numparticle),dxco(3,numparticle,ntime))
      allocate(vcor(3,numparticle),dvco(3,numparticle,ntime))
      !
      psize=msize(particle)
      !
      ! duplicat the particle array to local array
      !
      npart=0
      do jpart=1,psize
        !
        pa=>particle(jpart)
        !
        if(pa%new) cycle
        !
        npart=npart+1
        !
        xcor(:,npart)=pa%x(:)
        !
        dxco(:,npart,1)=pa%v(:)
        dxco(:,npart,2:ntime)=pa%dx(:,2:ntime)
        !
        vcor(:,npart)=pa%v(:)
        !
        dvco(:,npart,1)=pa%f(:)
        dvco(:,npart,2:ntime)=pa%dv(:,2:ntime)
        !
        if(npart==numparticle) exit
        !
      enddo
      !
      if (itimescheme.eq.1) then
         !>>> Euler
         vcor=gdt(itr)*tsro*dvco(:,:,1)+vcor
         !
         xcor=gdt(itr)*tsro*dxco(:,:,1)+xcor
         !
      elseif(itimescheme.eq.2) then
         !>>> Adam-Bashforth second order (AB2)
         !
         if(itime.eq.1 .and. irestart.eq.0) then
           ! Do first time step with Euler
           vcor=gdt(itr)*tsro*dvco(:,:,1)+vcor
           !
           xcor=gdt(itr)*tsro*dxco(:,:,1)+xcor
         else
           vcor=adt(itr)*dvco(:,:,1)+bdt(itr)*dvco(:,:,2)+vcor
           !
           xcor=adt(itr)*dxco(:,:,1)+bdt(itr)*dxco(:,:,2)+xcor
         endif
         dvco(:,:,2)=dvco(:,:,1)
         dxco(:,:,2)=dxco(:,:,1)
         !
      elseif(itimescheme.eq.3) then
         !>>> Adams-Bashforth third order (AB3)
         !
         ! Do first time step with Euler
         if(itime.eq.1.and.irestart.eq.0) then
            vcor=dt*tsro*dvco(:,:,1)+vcor
            !
            xcor=dt*tsro*dxco(:,:,1)+xcor
         elseif(itime.eq.2.and.irestart.eq.0) then
            ! Do second time step with AB2
            vcor=onepfive*dt*tsro*dvco(:,:,1)-half*dt*tsro*dvco(:,:,2)+vcor
            dvco(:,:,3)=dvco(:,:,2)
            !
            xcor=onepfive*dt*tsro*dxco(:,:,1)-half*dt*tsro*dxco(:,:,2)+xcor
            dxco(:,:,3)=dxco(:,:,2)
         else
            ! Finally using AB3
            vcor=adt(itr)*tsro*dvco(:,:,1)+bdt(itr)*tsro*dvco(:,:,2)+cdt(itr)*tsro*dvco(:,:,3)+vcor
            dvco(:,:,3)=dvco(:,:,2)
            !
            xcor=adt(itr)*tsro*dxco(:,:,1)+bdt(itr)*tsro*dxco(:,:,2)+cdt(itr)*tsro*dxco(:,:,3)+xcor
            dxco(:,:,3)=dxco(:,:,2)
         endif
         dvco(:,:,2)=dvco(:,:,1)
         dxco(:,:,2)=dxco(:,:,1)
         !
      elseif(itimescheme.eq.5) then
         !>>> Runge-Kutta (low storage) RK3
         if(itr.eq.1) then
            vcor=gdt(itr)*tsro*dvco(:,:,1)+vcor
            xcor=gdt(itr)*tsro*dxco(:,:,1)+xcor
         else
            vcor=adt(itr)*tsro*dvco(:,:,1)+bdt(itr)*tsro*dvco(:,:,2)+vcor
            xcor=adt(itr)*tsro*dxco(:,:,1)+bdt(itr)*tsro*dxco(:,:,2)+xcor
         endif
         dvco(:,:,2)=dvco(:,:,1)
         dxco(:,:,2)=dxco(:,:,1)
         !
      endif
      ! !
      ! put back from local array to particle array
      npart=0
      do jpart=1,psize
        !
        pa=>particle(jpart)
        !
        if(pa%new) cycle
        !
        npart=npart+1
        !
        if(isnan(vcor(1,npart)) .or. isnan(vcor(2,npart)) .or. isnan(vcor(3,npart))) then
          print*, '!! particle velocity NaN @ intt_particel'
          print*, pa%v(:),pa%f(:)
          stop
        endif
        !
        pa%v(:)=vcor(:,npart)
        pa%x(:)=xcor(:,npart)
        ! pa%x(3)=0.d0 ! 2d x-y computation
        !
        pa%dv(:,2:ntime)=dvco(:,npart,2:ntime)
        pa%dx(:,2:ntime)=dxco(:,npart,2:ntime)
        !
        if(npart==numparticle) exit
        !
      enddo
      !
      deallocate(xcor,dxco,vcor,dvco)
      !
      ! call partical_domain_check('bc_tgv')
      call partical_domain_check('out_disappear')
      !
      call partical_swap
      !
      part_time=part_time+ptime()-timebeg
      !
      total_num_part=psum(numparticle)
      !
      if(nrank==0) then
        if(total_num_part.ne.old_num_part) then
          print*,' ** number of particles changes from ',old_num_part,'->',total_num_part
        endif
        old_num_part=total_num_part
      endif
      !
      ! if(nrank==0) print*,particletime,time1
    enddo
    !
    ux0=ux1
    uy0=uy1
    uz0=uz1
    !
    time0=time1
    ! print*,nrank,'| number of particles:',numparticle
    ! !
    ! call mpistop
    !
    ! if(nrank==0) then
    !   print*,xpa(1),ypa(1),zpa(1),'|',ux_pa(1)
    ! endif
    !
  end subroutine intt_particel
  !+-------------------------------------------------------------------+
  ! The end of the subroutine intt_particel                            |
  !+-------------------------------------------------------------------+
  !
  function linintp(xx1,xx2,yy1,yy2,xx) result(yy)
    !
    real(8),intent(in) :: xx1,xx2,xx
    real(8),intent(in) ::  yy1(:,:,:),yy2(:,:,:)
    real(8) :: yy(1:size(yy1,1),1:size(yy1,2),1:size(yy1,3))
    !
    real(8) :: var1
    !
    var1=(xx-xx1)/(xx2-xx1)
    yy=(yy2-yy1)*var1+yy1
    !
    return
    !
  end function linintp
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to check if the particle is out of domain      |
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine partical_domain_check(mode)
    !
    use param,     only : nclx,ncly,nclz,xlx,yly,zlz
    !
    character(len=*),intent(in) :: mode
    !
    ! local data 
    integer :: jpart,npart,psize,jrank,npcanc,npcanc_totl
    type(partype),pointer :: pa
    real(mytype) :: timebeg
    !
    timebeg=ptime()
    !
    psize=msize(particle)
    !
    if(mode=='out_disappear') then
      !
      npart=0
      npcanc=0
      do jpart=1,psize
        !
        pa=>particle(jpart)
        !
        if(pa%new) cycle
        !
        npart=npart+1
        if(npart>numparticle) exit
        ! if the particle is out of domain, mark it and subscribe the 
        ! total number of particles
        !
        if(nclx .and. (pa%x(1)>xlx .or. pa%x(1)<0)) then
          call pa%reset()
          npcanc=npcanc+1
          cycle
        endif
        !
        if(ncly .and. (pa%x(2)>yly .or. pa%x(2)<0)) then
          call pa%reset()
          npcanc=npcanc+1
          cycle
        endif
        !
        if(nclz .and. (pa%x(3)>zlz .or. pa%x(3)<0)) then
          call pa%reset()
          npcanc=npcanc+1
          cycle
        endif
        !
        if( pa%x(1)>=lxmin(nrank) .and. pa%x(1)<lxmax(nrank) .and. &
            pa%x(2)>=lymin(nrank) .and. pa%x(2)<lymax(nrank) .and. &
            pa%x(3)>=lzmin(nrank) .and. pa%x(3)<lzmax(nrank) ) then
          continue
        else
          !
          pa%swap=.true.
          !
          do jrank=0,nproc-1
            !
            ! to find which rank the particle are moving to and 
            ! mark
            if(jrank==nrank) cycle
            !
            if( pa%x(1)>=lxmin(jrank) .and. pa%x(1)<lxmax(jrank) .and. &
                pa%x(2)>=lymin(jrank) .and. pa%x(2)<lymax(jrank) .and. &
                pa%x(3)>=lzmin(jrank) .and. pa%x(3)<lzmax(jrank) ) then
              !
              pa%rank2go=jrank
              !
              exit
              !
            endif
            !
          enddo
          !
        endif
        !
      enddo
      !
      numparticle=numparticle-npcanc
      !
      npcanc_totl=psum(npcanc)
      if(nrank==0 .and. npcanc_totl>0) print*,' ** ',npcanc_totl,        &
                                     ' particles are moving out of domain'
      !
    elseif(mode=='bc_channel') then
      !
      npart=0
      npcanc=0
      do jpart=1,psize
        !
        pa=>particle(jpart)
        !
        if(pa%new) cycle
        !
        npart=npart+1
        !
        if(npart>numparticle) exit
        ! if the particle is out of domain, mark it and subscribe the 
        ! total number of particles
        !
        if(nclx) then
          !
          if(pa%x(1)>xlx) then
            pa%x(1)=pa%x(1)-xlx
          endif
          !
          if(pa%x(1)<0) then
            pa%x(1)=pa%x(1)+xlx
          endif
          !
        endif
        !
        if(ncly) then
          !
          if(pa%x(2)>yly) then
            pa%x(2)=pa%x(2)-yly
          endif 
          !
          if(pa%x(2)<0) then
            pa%x(2)=pa%x(2)+yly
          endif
          !
        else
          !
          ! reflect particles back in the domain.
          if(pa%x(2)>yly) then
            pa%x(2)=2.d0*yly-pa%x(2)
          endif
          if(pa%x(2)<0) then
            pa%x(2)=-pa%x(2)
          endif
          !
        endif
        !
        if(nclz) then
          !
          if(pa%x(3)>zlz) then
            pa%x(3)=pa%x(3)-zlz
          endif 
          !
          if(pa%x(3)<0) then
            pa%x(3)=pa%x(3)+zlz
          endif
          !
        endif
        !
        if(pa%x(1)>xlx .or. pa%x(1)<0 .or. &
           pa%x(2)>yly .or. pa%x(2)<0 .or. &
           pa%x(3)>zlz .or. pa%x(3)<0 ) then
            print*,' !! waring, the particle still moves out of domain 1'
            print*,nrank,jpart,'x:',pa%x(:),'v:',pa%v(:),'vf:',pa%vf(:),'f:',pa%f(:)
            stop
        endif
        !
        if( pa%x(1)>=lxmin(nrank) .and. pa%x(1)<lxmax(nrank) .and. &
            pa%x(2)>=lymin(nrank) .and. pa%x(2)<lymax(nrank) .and. &
            pa%x(3)>=lzmin(nrank) .and. pa%x(3)<lzmax(nrank) ) then
          continue
        else
          !
          pa%swap=.true.
          !
          do jrank=0,nproc-1
            !
            ! to find which rank the particle are moving to and 
            ! mark
            if(jrank==nrank) cycle
            !
            if( pa%x(1)>=lxmin(jrank) .and. pa%x(1)<lxmax(jrank) .and. &
                pa%x(2)>=lymin(jrank) .and. pa%x(2)<lymax(jrank) .and. &
                pa%x(3)>=lzmin(jrank) .and. pa%x(3)<lzmax(jrank) ) then
              !
              pa%rank2go=jrank
              !
              exit
              !
            endif
            !
          enddo
          !
        endif
        !
      enddo
      !
    elseif(mode=='bc_tgv') then
      !
      npart=0
      npcanc=0
      do jpart=1,psize
        !
        pa=>particle(jpart)
        !
        if(pa%new) cycle
        !
        npart=npart+1
        !
        if(npart>numparticle) exit
        ! if the particle is out of domain, mark it and subscribe the 
        ! total number of particles
        !
        if(nclx) then
          !
          if(pa%x(1)>xlx) then
            pa%x(1)=pa%x(1)-xlx
          endif
          !
          if(pa%x(1)<0) then
            pa%x(1)=pa%x(1)+xlx
          endif
          !
        else
          stop ' 1 @ partical_domain_check'
        endif
        !
        if(ncly) then
          !
          if(pa%x(2)>yly) then
            pa%x(2)=pa%x(2)-yly
          endif 
          !
          if(pa%x(2)<0) then
            pa%x(2)=pa%x(2)+yly
          endif
          !
        else
          print*,' ** particles x: ',pa%x(:)
          stop ' 2 @ partical_domain_check'
        endif
        !
        if(nclz) then
          !
          if(pa%x(3)>zlz) then
            pa%x(3)=pa%x(3)-zlz
          endif 
          !
          if(pa%x(3)<0) then
            pa%x(3)=pa%x(3)+zlz
          endif
          !
        else
          stop ' 3 @ partical_domain_check'
        endif
        !
        if(pa%x(1)>xlx .or. pa%x(1)<0 .or. &
           pa%x(2)>yly .or. pa%x(2)<0 .or. &
           pa%x(3)>zlz .or. pa%x(3)<0 ) then
            print*,' !! waring, the particle still moves out of domain 2'
            print*,nrank,jpart,'x:',pa%x(:),'v:',pa%v(:),'vf:',pa%vf(:),'f:',pa%f(:)
            stop
        endif
        !
        if( pa%x(1)>=lxmin(nrank) .and. pa%x(1)<lxmax(nrank) .and. &
            pa%x(2)>=lymin(nrank) .and. pa%x(2)<lymax(nrank) .and. &
            pa%x(3)>=lzmin(nrank) .and. pa%x(3)<lzmax(nrank) ) then
          continue
        else
          !
          pa%swap=.true.
          !
          do jrank=0,nproc-1
            !
            ! to find which rank the particle are moving to and 
            ! mark
            if(jrank==nrank) cycle
            !
            if( pa%x(1)>=lxmin(jrank) .and. pa%x(1)<lxmax(jrank) .and. &
                pa%x(2)>=lymin(jrank) .and. pa%x(2)<lymax(jrank) .and. &
                pa%x(3)>=lzmin(jrank) .and. pa%x(3)<lzmax(jrank) ) then
              !
              pa%rank2go=jrank
              !
              exit
              !
            endif
            !
          enddo
          !
        endif
        !
      enddo
      !
    else
      !
      stop ' !! mode not defined @ partical_domain_check!!'
      !
    endif
    !
    part_dmck_time=part_dmck_time+ptime()-timebeg
    ! print*,nrank,'|',numparticle
    ! do jpart=1,psize
    !   !
    !   pa=>particle(jpart)
    !   !
    !   if(pa%swap) then
    !     !
    !     write(*,'(3(A,1X,I0))')' ** particle',jpart,' moves from rank', &
    !                                  pa%rankinn,' to rank',pa%rank2go
    !     !
    !   endif
    !   !
    ! enddo
    !
  end subroutine partical_domain_check
  !+-------------------------------------------------------------------+
  ! The end of the subroutine partical_domain_check                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to swap particle infomation between ranks      |
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine partical_swap
    !
    use param,     only : nclx,ncly,nclz,xlx,yly,zlz
    !
    ! local data 
    integer :: p,psize,jrank,jpart,npart,n,newsize
    type(partype),pointer :: pa
    integer :: nsend(0:nproc-1),nrecv(0:nproc-1),nsend_total
    !+------------+--------------------------------+
    !| nsendtable | the table to record how much   |
    !|            | particle to send to a rank     |
    !+------------+--------------------------------+
    integer :: pr(0:nproc-1,1:numparticle)
    type(partype),allocatable :: pa2send(:),pa2recv(:)
    real(mytype) :: timebeg,tvar1,tvar11,tvar2,tvar3,tvar4
    !
    timebeg=ptime()
    !
    psize=msize(particle)
    !
    nsend=0
    !
    n=0
    pr=0
    npart=0
    !
    do jpart=1,psize
      !
      pa=>particle(jpart)
      !
      if(pa%new) cycle
      !
      ! to find out how many particle to send to which ranks
      if(pa%swap) then
        !
        n=n+1
        !
        nsend(pa%rank2go)=nsend(pa%rank2go)+1
        !
        pr(pa%rank2go,nsend(pa%rank2go))=jpart
        !
      endif
      !
      npart=npart+1
      !
      if(npart==numparticle) exit
      !
    enddo
    !
    nsend_total=n
    !
    tvar1=ptime()
    count_time=count_time+tvar1-timebeg
    !
    ! do jrank=0,nproc-1
    !   if(nsend(jrank)>0) then
    !     print*,' **',nsend(jrank),'particles is moving ',nrank,'->',jrank
    !   endif
    ! enddo
    !
    ! synchronize recv table according to send table
    nrecv=ptabupd(nsend)
    !
    tvar11=ptime()
    !
    table_share_time=table_share_time+tvar11-tvar1
    !
    !
    ! to establish the buffer of storing particels about to send
    if(nsend_total>0) then
      !
      allocate(pa2send(1:nsend_total))
      !
      n=0
      do jrank=0,nproc-1
        !
        do jpart=1,nsend(jrank)
          !
          n=n+1
          !
          p=pr(jrank,jpart)
          !
          pa2send(n)=particle(p)
          !
          call particle(p)%reset()
          !
          numparticle=numparticle-1
          !
        enddo
        !
      enddo
      !
    endif 
    !
    tvar2=ptime()
    data_pack_time=data_pack_time+tvar2-tvar11
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! swap particle among ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call pa2a(pa2send,pa2recv,nsend,nrecv)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of swap particle among ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    tvar3=ptime()
    mpi_comm_time=mpi_comm_time+tvar3-tvar2
    !
    call particle_add(particle,pa2recv,n)
    !
    ! do jrank=0,nproc-1
    !   if(n>0) then
    !     print*,nrank,'| add ',n,'particles'
    !   endif
    ! enddo
    ! now add the received particle in to the array, dynamically
    ! if(numparticle+msize(pa2recv)>psize) then
    !   !
    !   ! expand the particle array
    !   newsize=max(numparticle+msize(pa2recv),numparticle+100)
    !   !
    !   call mextend(particle,newsize)
    !   !
    ! endif
    ! !
    ! n=0
    ! do jpart=1,msize(particle)
    !   !
    !   pa=>particle(jpart)
    !   !
    !   ! the particle is free for re-assigning
    !   if(pa%new) then
    !     !
    !     if(n>=msize(pa2recv)) exit
    !     !
    !     n=n+1
    !     !
    !     pa=pa2recv(n)
    !     pa%new=.false.
    !     !
    !   endif
    !   !
    ! enddo
    ! !
    numparticle=numparticle+n
    !
    tvar4=ptime()
    data_unpack_time=data_unpack_time+tvar4-tvar3
    !
    part_comm_time=part_comm_time+ptime()-timebeg
    !
  end subroutine partical_swap
  !+-------------------------------------------------------------------+
  ! The end of the subroutine partical_swap                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to add particles to the current particle arrary|
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine particle_add(particle_cur,particle_new,num_part_incr)
    !
    ! arguments
    type(partype),intent(inout),allocatable,target :: particle_cur(:)
    type(partype),intent(in),allocatable :: particle_new(:)
    integer,intent(out) :: num_part_incr
    !
    ! local data
    integer :: psize,newsize,n,jpart
    type(partype),pointer :: pa
    !
    psize=msize(particle_cur)
    !
    ! now add the received particle in to the array, dynamically
    if(numparticle+msize(particle_new)>psize) then
      !
      ! expand the particle array
      newsize=max(numparticle+msize(particle_new),numparticle+100)
      !
      call mextend(particle_cur,newsize)
      !
      psize=newsize
    endif
    !
    n=0
    do jpart=1,psize
      !
      ! print*,nrank,'|',jpart
      pa=>particle_cur(jpart)
      !
      ! the particle is free for re-assigning
      if(pa%new) then
        !
        if(n>=msize(particle_new)) exit
        !
        n=n+1
        !
        pa=particle_new(n)
        pa%new=.false.
        !
      endif
      !
    enddo
    !
    ! print*,nrank,'|',n,newsize
    num_part_incr=n
    !
  end subroutine particle_add
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_add                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to write particle coordinate.                  |
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine write_particle
    !
    use param, only : t
    !
    ! local data
    integer :: p,psize,total_num_part
    logical,save :: firstcal=.true.
    !
    if(firstcal .and. numparticle>0) then
      !
      write(rankname,'(i4.4)')nrank
      !
      open(18,file='./data/particle'//rankname//'.dat')
      write(18,'(A)')'VARIABLES = "x" "y" "z" "nrank" '
      close(18)
      print*,' create particle',rankname,'.dat'
      !
      firstcal=.false.
      !
    endif
    !
    psize=msize(particle)
    !
    if(numparticle>0) then
      open(18,file='./data/particle'//rankname//'.dat',position="append")
      write(18,'(A)')'ZONE T="ZONE 001"'
      write(18,'(A,I0,A)')'I=',numparticle,', J=1, K=1, ZONETYPE=Ordered'
      write(18,'(A,E13.6E2)')'STRANDID=1, SOLUTIONTIME=',t
      write(18,'(A)')'DATAPACKING=POINT'
      do p=1,psize
        if(particle(p)%new) cycle
        write(18,'(3(1X,E15.7E3),1X,I0)')particle(p)%x(1),particle(p)%x(2),particle(p)%x(3),particle(p)%rankinn
      enddo
      close(18)
      ! print*,' << ./data/particle',rankname,'.dat'
    endif
    !
    total_num_part=psum(numparticle)
    !
    if(nrank==0) print*,' ** total number of particles is:',total_num_part
    ! print*,nrank,'| number of particles:',numparticle
    !
  end subroutine write_particle
  !+-------------------------------------------------------------------+
  ! The end of the subroutine write_particle                           |
  !+-------------------------------------------------------------------+
  !
  subroutine h5read_particle(particle_new,particle_size)
    !
    use param, only : t,itime
    !
    ! argument
    type(partype),intent(out),allocatable :: particle_new(:)
    integer,intent(out) :: particle_size
    !
    ! local data
    character(len=5) :: num
    character(len=32) :: file2read
    integer :: nstep,psize,jp
    real(mytype) :: time
    real(mytype),allocatable :: xpart(:),ypart(:),zpart(:)
    !
    if(nrank==0) then
      !
      write(num,'(I5.5)') particle_file_numb
      !
      file2read='./data/particle'//num//'.h5'
      !
      ! call h5sread(var=nstep,varname='itime',filename=trim(file2read),explicit=.true.)
      ! !
      ! if(nstep .ne. itime) then
      !   print*,' the itime from the particle file ',trim(file2read),   &
      !                             'not consistent with the restart file'
      !   stop ' !! ERROR @ h5read_particle'
      ! endif
      ! !
      ! call h5sread(var=time,varname='time',filename=trim(file2read),explicit=.true.)
      ! print*,' ** time: ',time
      !
      psize=h5getdim3d(varname='x',filenma=trim(file2read))
      print*,' ** number of particles: ',psize
      !
      allocate(xpart(psize),ypart(psize),zpart(psize))
      call h5sread(var=xpart,varname='x',dim=psize,filename=trim(file2read),explicit=.true.)
      call h5sread(var=ypart,varname='y',dim=psize,filename=trim(file2read),explicit=.true.)
      call h5sread(var=zpart,varname='z',dim=psize,filename=trim(file2read),explicit=.true.)
      !
      allocate(particle_new(1:psize))
      !
      do jp=1,psize
        !
        call particle_new(jp)%init()
        !
        particle_new(jp)%x(1)=xpart(jp)
        particle_new(jp)%x(2)=ypart(jp)
        particle_new(jp)%x(3)=zpart(jp)
        !
        particle_new(jp)%new=.false.
        !
      enddo
      !
      particle_size=psize
      !
    else
      particle_size=0
    endif
    !
  end subroutine h5read_particle
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to write particles via HDF5.                   |
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-09-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine h5write_particle()
    !
    use param, only : t,itime,irestart
    !
    ! local data
    character(len=5) :: num
    real(mytype),allocatable :: xpart(:),ypart(:),zpart(:),                 &
                           upart(:),vpart(:),wpart(:)
    integer :: psize,total_num_part,p,j
    integer :: rank2coll
    character(len=32) :: file2write
    real(mytype) :: timebeg
    logical :: fexists
    !
    logical,save :: init=.true.
    !
    timebeg=ptime()
    !
    if(init) then
      !
      ! write head of the xdmf file
      if(nrank==0) then
        !
        inquire(file='./visu_particle.xdmf', exist=fexists)
        !
        if(fexists .and. irestart.ne.0) then
          continue
        else
          open(22,file='./visu_particle.xdmf',form='formatted')
          write(22,'(A)')'<?xml version="1.0" encoding="UTF-8" ?>'
          write(22,'(A)')'<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
          write(22,'(A)')'<Domain>'
          write(22,'(A)')'  <Grid Name="Particle" GridType="Collection" CollectionType="Temporal">'
          write(22,'(A)')'  </Grid>'
          write(22,'(A)')'</Domain>'
          write(22,'(A)')'</Xdmf>'
          close(22)
          print*,' << visu_particle.xdmf'
        endif
        !
      endif
      !
      init=.false.
      !
    endif
    !
    if(numparticle>0) then
      !
      allocate(xpart(numparticle),ypart(numparticle),zpart(numparticle))
      allocate(upart(numparticle),vpart(numparticle),wpart(numparticle))
      !
      psize=msize(particle)
      j=0
      do p=1,psize
        !
        if(particle(p)%new) cycle
        !
        j=j+1
        xpart(j)=particle(p)%x(1)
        ypart(j)=particle(p)%x(2)
        zpart(j)=particle(p)%x(3)
        !
        upart(j)=particle(p)%v(1)
        vpart(j)=particle(p)%v(2)
        wpart(j)=particle(p)%v(3)
        !
      enddo
      !
      rank2coll=nrank
      !
    else
      rank2coll=-1
    endif
    !
    total_num_part=psum(numparticle)
    !
    call subcomm_group(rank2coll,mpi_comm_particle,mpi_rank_part,mpi_size_part)
    !
    ! print*,' ** size of the particel comm:',mpi_size_part
    !
    particle_file_numb=particle_file_numb+1
    !
    if(rank2coll>=0) then
      !
      write(num,'(I5.5)') particle_file_numb
      !
      file2write='./data/particle'//num//'.h5'
      call h5io_init(filename=trim(file2write),mode='writ',            &
                                                 comm=mpi_comm_particle)
      !
      call h5write(varname='itime',var=itime)
      call h5write(varname='time',var=t)
      call h5write(varname='x',var=xpart,total_size=total_num_part,    &
                                               comm=mpi_comm_particle, &
                                          comm_size=mpi_size_part,     &
                                          comm_rank=mpi_rank_part)
      call h5write(varname='y',var=ypart,total_size=total_num_part,    &
                                               comm=mpi_comm_particle, &
                                          comm_size=mpi_size_part,     &
                                          comm_rank=mpi_rank_part)
      call h5write(varname='z',var=zpart,total_size=total_num_part,    &
                                               comm=mpi_comm_particle, &
                                          comm_size=mpi_size_part,     &
                                          comm_rank=mpi_rank_part)
      call h5write(varname='u',var=upart,total_size=total_num_part,    &
                                               comm=mpi_comm_particle, &
                                          comm_size=mpi_size_part,     &
                                          comm_rank=mpi_rank_part)
      call h5write(varname='v',var=vpart,total_size=total_num_part,    &
                                               comm=mpi_comm_particle, &
                                          comm_size=mpi_size_part,     &
                                          comm_rank=mpi_rank_part)
      call h5write(varname='w',var=wpart,total_size=total_num_part,    &
                                               comm=mpi_comm_particle, &
                                          comm_size=mpi_size_part,     &
                                          comm_rank=mpi_rank_part)
      !
      deallocate(xpart,ypart,zpart,upart,vpart,wpart)
      !
      call h5io_end
      !
      if(mpi_rank_part==0) then
        !
        print*,' << ',trim(file2write)
        !
        open(22,file='./visu_particle.xdmf',form='formatted',position="append")
        backspace(22)
        backspace(22)
        backspace(22)
        write(22,'(A,A,A)')    '    <Grid Name=" particle',num,'" GridType="Uniform">'
        write(22,'(A,F12.6,A)')'    <Time Value="',t,'" />'
        write(22,'(A,I0,A)')   '    <Topology Name="ParticleTopo" TopologyType="PolyVertex" NumberOfElements="',total_num_part,'"/>'
        write(22,'(A)')        '    <Geometry GeometryType="X_Y_Z">'
        write(22,'(A,I0,A)')   '    <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',total_num_part,'">'
        write(22,'(A,A,A)')    '       ',trim(file2write),':x'
        write(22,'(A)')        '    </DataItem>'
        write(22,'(A,I0,A)')   '    <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',total_num_part,'">'
        write(22,'(A,A,A)')    '       ',trim(file2write),':y'
        write(22,'(A)')        '    </DataItem>'
        write(22,'(A,I0,A)')   '    <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',total_num_part,'">'
        write(22,'(A,A,A)')    '       ',trim(file2write),':z'
        write(22,'(A)')        '    </DataItem>'
        write(22,'(A)')        '    </Geometry>'
        write(22,'(A)')        '    </Grid>' 
        write(22,'(A)')        '  </Grid>'
        write(22,'(A)')        '</Domain>'
        write(22,'(A)')        '</Xdmf>'
        !
        close(22)
        print*,' << visu_particle.xdmf'
      endif
      !
    endif
    !
    h5io_time=h5io_time+ptime()-timebeg
    !
    if(nrank==0) print*,' ** total number of particles is:',total_num_part
    !
  end subroutine h5write_particle
  !
  subroutine pswap_yz(varin,varout)
    !
    use variables, only : p_row, p_col
    use decomp_2d, only : xsize
    use param,     only : nclx,ncly,nclz
    use mpi
    !
    ! arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: varin
    real(mytype),intent(out) :: varout(1:xsize(1)+1,0:xsize(2)+1,0:xsize(3)+1)
    !
    ! local data
    integer :: i,j,k,n,jrk,krk,ncou,mpitag,ierr
    integer :: status(mpi_status_size) 
    integer,allocatable,save :: mrank(:,:)
    integer,save :: upper,lower,front,bback
    logical,save :: init=.true.
    real(mytype),dimension(:,:),allocatable :: sbuf1,sbuf2,rbuf1,rbuf2
    !
    if(init) then
      !
      allocate(mrank(p_row,p_col))
      !
      n=-1
      do j=1,p_row
      do k=1,p_col
        !
        n=n+1
        mrank(j,k)=n
        !
        if(nrank==n) then
          jrk=j
          krk=k
        endif
        !
      enddo
      enddo
      !
      if(jrk==1) then
        upper=mrank(jrk+1,krk)
        !
        if(ncly) then
          lower=mrank(p_row,krk)
        else
          lower=mpi_proc_null
        endif
        !
      elseif(jrk==p_row) then
        !
        if(ncly) then
          upper=mrank(1,krk)
        else
          upper=mpi_proc_null
        endif
        !
        lower=mrank(jrk-1,krk)
      else
        upper=mrank(jrk+1,krk)
        lower=mrank(jrk-1,krk)
      endif
      !
      if(p_col==1) then
        front=mpi_proc_null
        bback=mpi_proc_null
      else
        !
        if(krk==1) then
          front=mrank(jrk,krk+1)
          !
          if(nclz) then
            bback=mrank(jrk,p_col)
          else
            bback=mpi_proc_null
          endif
          !
        elseif(krk==p_col) then
          if(nclz) then
            front=mrank(jrk,1)
          else
            front=mpi_proc_null
          endif
          !
          bback=mrank(jrk,krk-1)
        else
          front=mrank(jrk,krk+1)
          bback=mrank(jrk,krk-1)
        endif
      endif
      !
      ! print*,p_row,p_col
      ! print*,nrank,'-',jrk,krk,':',upper,lower,front,bback,mpi_proc_null
      !
      init=.false.
      !
    endif
    !
    varout(1:xsize(1),1:xsize(2),1:xsize(3))=varin(1:xsize(1),1:xsize(2),1:xsize(3))
    !
    mpitag=1000
    !
    ! send & recv in the z direction
    !
    if(xsize(3)==1) then
      varout(1:xsize(1),1:xsize(2),0)=varout(1:xsize(1),1:xsize(2),1)
      varout(1:xsize(1),1:xsize(2),2)=varout(1:xsize(1),1:xsize(2),1)
    else
      !
      allocate(sbuf1(1:xsize(1),1:xsize(2)),sbuf2(1:xsize(1),1:xsize(2)),&
               rbuf1(1:xsize(1),1:xsize(2)),rbuf2(1:xsize(1),1:xsize(2)))
      !
      if(front .ne. mpi_proc_null) then
        sbuf1(1:xsize(1),1:xsize(2))=varout(1:xsize(1),1:xsize(2),xsize(3))
      endif
      if(bback .ne. mpi_proc_null) then
        sbuf2(1:xsize(1),1:xsize(2))=varout(1:xsize(1),1:xsize(2),1)
      endif
      !
      ncou=xsize(1)*xsize(2)
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,front, mpitag,             &
                        rbuf1,ncou,mpi_real8,bback, mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,bback, mpitag,            &
                        rbuf2,ncou,mpi_real8,front, mpitag,            &
                                               mpi_comm_world,status,ierr)
      !
      if(bback .ne. mpi_proc_null) then
        varout(1:xsize(1),1:xsize(2),0)=rbuf1(1:xsize(1),1:xsize(2))
      endif
      if(front .ne. mpi_proc_null) then
        varout(1:xsize(1),1:xsize(2),xsize(3)+1)=rbuf2(1:xsize(1),1:xsize(2))
      endif
      !
      deallocate(sbuf1,sbuf2,rbuf1,rbuf2)
      !
    endif
    !
    ! end of Message passing in the z direction
    !
    ! 
    ! send & recv in the y direction
    allocate( sbuf1(1:xsize(1),0:xsize(3)+1),                          &
              sbuf2(1:xsize(1),0:xsize(3)+1),                          &
              rbuf1(1:xsize(1),0:xsize(3)+1),                          &
              rbuf2(1:xsize(1),0:xsize(3)+1) )

    if(upper .ne. mpi_proc_null) then
      sbuf1(1:xsize(1),0:xsize(3)+1)=varout(1:xsize(1),xsize(2),0:xsize(3)+1)
    endif
    if(lower .ne. mpi_proc_null) then
      sbuf2(1:xsize(1),0:xsize(3)+1)=varout(1:xsize(1),1,0:xsize(3)+1)
    endif
    !
    ncou=xsize(1)*(xsize(3)+2)
    !
    ! Message passing
    call mpi_sendrecv(sbuf1,ncou,mpi_real8,upper, mpitag,             &
                      rbuf1,ncou,mpi_real8,lower, mpitag,             &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(sbuf2,ncou,mpi_real8,lower, mpitag,            &
                      rbuf2,ncou,mpi_real8,upper, mpitag,            &
                                             mpi_comm_world,status,ierr)
    !
    if(upper .ne. mpi_proc_null) then
      varout(1:xsize(1),xsize(2)+1,0:xsize(3)+1)=rbuf2(1:xsize(1),0:xsize(3)+1)
    endif
    if(lower .ne. mpi_proc_null) then
      varout(1:xsize(1),0,0:xsize(3)+1)=rbuf1(1:xsize(1),0:xsize(3)+1)
    endif
    !
    deallocate(sbuf1,sbuf2,rbuf1,rbuf2)
    ! send & recv in the y direction
    !
    if(nclx) then
      varout(xsize(1)+1,:,:)=varout(1,:,:)
    endif
    !
    return
    !
  end subroutine pswap_yz
  !
  !+-------------------------------------------------------------------+
  !| this subroutine clean superfluous elements in a array             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-Nov-2018  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine mclean_mytype(var,n)
    !
    ! arguments
    real(mytype),allocatable,intent(inout) :: var(:)
    integer,intent(in) :: n
    !
    ! local data
    real(mytype),allocatable :: buffer(:)
    integer :: m
    logical :: lefc
    !
    if(.not.allocated(var)) return
    !
    if(n<=0) then
      deallocate(var)
      return
    endif
    !
    ! clean
    allocate(buffer(n))
    !
    buffer(1:n)=var(1:n)
    !
    deallocate(var)
    !
    call move_alloc(buffer,var)
    !
  end subroutine mclean_mytype
  !!
  subroutine mclean_particle(var,n)
    !
    ! arguments
    type(partype),allocatable,intent(inout) :: var(:)
    integer,intent(in) :: n
    !
    ! local data
    type(partype),allocatable :: buffer(:)
    integer :: m
    logical :: lefc
    !
    if(.not.allocated(var)) return
    !
    if(n<=0) then
      deallocate(var)
      return
    endif
    !
    ! clean
    allocate(buffer(n))
    !
    buffer(1:n)=var(1:n)
    !
    deallocate(var)
    !
    call move_alloc(buffer,var)
    !
  end subroutine mclean_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mclean.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this function is to retune the size of a array                    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  pure function size_particle(var) result(nsize)
    !
    type(partype),allocatable,intent(in) :: var(:)
    integer :: nsize
    !
    if(allocated(var)) then
      nsize=size(var)
    else
      nsize=0
    endif
    !
    return
    !
  end function size_particle
  !
  pure function size_integer(var) result(nsize)
    !
    integer,allocatable,intent(in) :: var(:)
    integer :: nsize
    !
    if(allocated(var)) then
      nsize=size(var)
    else
      nsize=0
    endif
    !
    return
    !
  end function size_integer
  !+-------------------------------------------------------------------+
  !| The end of the subroutine size_particle.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is to swap particles via alltoall.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine pa2a_particle(datasend,datarecv,sendtabl,recvtabl)
    !
    use mpi
    !
    ! arguments
    type(partype),allocatable,intent(in) ::  datasend(:)
    type(partype),allocatable,intent(out) :: datarecv(:)
    integer,intent(in) :: sendtabl(0:),recvtabl(0:)
    !
    ! local data
    integer :: ierr,recvsize,jrank,jpart,jc,nindsize
    integer,allocatable :: senddispls(:),recvdispls(:)
    real(mytype),allocatable :: r8send(:,:),r8resv(:,:)
    !
    integer,save :: newtype
    !
    logical,save :: firstcal=.true.
    !
    real(mytype) :: timebeg
    !
    timebeg=ptime()
    !
    if(firstcal) then
      call mpi_type_contiguous(6,mpi_real8,newtype,ierr)
      call mpi_type_commit(newtype,ierr)
      firstcal=.false.
    endif
    !
    allocate(senddispls(0:nproc-1),recvdispls(0:nproc-1))
    !
    senddispls=0
    recvdispls=0
    do jrank=1,nproc-1
      senddispls(jrank)=senddispls(jrank-1)+sendtabl(jrank-1)
      recvdispls(jrank)=recvdispls(jrank-1)+recvtabl(jrank-1)
    enddo
    recvsize=recvdispls(nproc-1)+recvtabl(nproc-1)
    !
    nindsize=msize(datasend)
    !
    allocate(r8send(6,nindsize))
    allocate(r8resv(6,recvsize))
    !
    r8resv=0.d0
    !
    do jpart=1,nindsize
      r8send(1,jpart)=datasend(jpart)%x(1)
      r8send(2,jpart)=datasend(jpart)%x(2)
      r8send(3,jpart)=datasend(jpart)%x(3)
      r8send(4,jpart)=datasend(jpart)%v(1)
      r8send(5,jpart)=datasend(jpart)%v(2)
      r8send(6,jpart)=datasend(jpart)%v(3)
    enddo
    !
    call mpi_alltoallv(r8send, sendtabl, senddispls, newtype, &
                       r8resv, recvtabl, recvdispls, newtype, &
                       mpi_comm_world, ierr)
    !
    allocate(datarecv(recvsize))
    !
    jc=0
    do jrank=0,nproc-1
      do jpart=1,recvtabl(jrank)
        !
        jc=jc+1
        !
        call datarecv(jc)%init()
        !
        datarecv(jc)%x(1)=r8resv(1,jc)
        datarecv(jc)%x(2)=r8resv(2,jc)
        datarecv(jc)%x(3)=r8resv(3,jc)
        datarecv(jc)%v(1)=r8resv(4,jc)
        datarecv(jc)%v(2)=r8resv(5,jc)
        datarecv(jc)%v(3)=r8resv(6,jc)
        !
        ! print*,nrank,'|',datarecv(jc)%x(1),datarecv(jc)%x(2),datarecv(jc)%x(3)
        !
      enddo
    enddo
    !
    a2a_time=a2a_time+ptime()-timebeg
    !
  end subroutine pa2a_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pa2a_particle.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is to expand an array.                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine extend_particle(var,n)
    !
    ! arguments
    type(partype),allocatable,intent(inout) :: var(:)
    integer,intent(in) :: n
    !
    ! local data
    type(partype),allocatable :: buffer(:)
    integer :: m,jpart
    !
    if(.not. allocated(var)) then
      allocate(var(n))
      m=0
    else
      !
      m=size(var)
      !
      call move_alloc(var, buffer)
      !
      allocate(var(n))
      var(1:m)=buffer(1:m)
      !
    endif
    !
    ! initilise newly allocated particles
    do jpart=m+1,n
      call var(jpart)%init()
    enddo
    !
    return
    !
  end subroutine extend_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pa2a_particle.                          |
  !+-------------------------------------------------------------------+
  !
end module partack
!+---------------------------------------------------------------------+
! The end of the module partack                                        |
!+---------------------------------------------------------------------+