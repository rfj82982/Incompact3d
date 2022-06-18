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
  !
  implicit none
  !
  interface psum
    module procedure psum_mytype_ary
    module procedure psum_integer
  end interface
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
  interface ptabupd
    module procedure ptable_update_int_arr
  end interface ptabupd
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
    real(mytype) :: x,y,z,u,v,w
    integer :: id,rankinn,rank2go
    logical :: swap,new
    !+------------------+------------------------------------------+
    !|            x,y,z | spatial coordinates of particle          |
    !|            u,v,w | velocity of particle                     |
    !|               id | the identification of particle           |
    !|          rankinn | the mpi rank which the particle is in    |
    !|          rank2go | the mpi rank which the particle will be  |
    !+------------------+------------------------------------------+
    !
    contains
    !
    procedure :: init => init_one_particle
    procedure :: reset => init_one_particle
    !
  end type partype
  !
  type(partype),allocatable,target :: particle(:)
  !
  logical :: lpartack
  integer :: numparticle,ipartiout,ipartiadd
  real(mytype) :: partirange(6)
  integer :: numpartix(3)
  real(mytype),allocatable,dimension(:) :: lxmin,lxmax,lymin,lymax,lzmin,lzmax
  real(mytype),allocatable,dimension(:) :: xpa,ypa,zpa
  real(mytype),allocatable,dimension(:) :: ux_pa,uy_pa,uz_pa
  character(len=4) :: rankname
  real(8) :: part_time,part_comm_time,part_vel_time,part_dmck_time,a2a_time
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
    class(partype),target :: pa
    !
    pa%swap=.false.
    pa%new =.true.
    !
    pa%rankinn=nrank
    pa%rank2go=nrank
    !
    pa%x=0.d0; pa%u=0.d0
    pa%y=0.d0; pa%v=0.d0
    pa%z=0.d0; pa%w=0.d0
    !
  end subroutine init_one_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine initmesg.                               |
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
      write(*,*) '      time particles com :',real(part_comm_time,4)
      write(*,*) '      time alltoall comm :',real(a2a_time,4)
      write(*,*) '      time domain search :',real(part_dmck_time,4)
    endif
    !
  end subroutine partcle_report
  !+-------------------------------------------------------------------+
  !| The end of the subroutine partcle_report.                         |
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
    use param,     only : xlx,yly,zlz
    !
    ! local data
    integer :: i,j,k,p
    integer :: big_num
    real(mytype) :: dx,dy,dz
    !
    !
    p=0
    !
    if(nrank==0) then
      !
      numparticle=numpartix(1)*numpartix(2)*numpartix(3)
      allocate(particle(1:numparticle))
      !
      do k=1,numpartix(3)
      do j=1,numpartix(2)
      do i=1,numpartix(1)
        !
        dy=(partirange(4)-partirange(3))/real(numpartix(2)-1,mytype)
        dz=(partirange(6)-partirange(5))/real(numpartix(3)-1,mytype)
        !
        p=p+1
        !
        call particle(p)%init()
        !
        particle(p)%x=0.1d0
        particle(p)%y=dy*real(j-1,mytype)+partirange(3)
        particle(p)%z=dz*real(k-1,mytype)+partirange(5)
        !
        particle(p)%new=.false.
        !
      enddo
      enddo
      enddo
      !
      ! numparticle=p
      ! !
      ! call mclean(particle,numparticle)
      !
    endif
    !
    ! p=0
    ! do k=1,numpartix(3)
    ! do j=1,numpartix(2)
    ! do i=1,numpartix(1)
    !   !
    !   dy=(lymax(nrank)-lymin(nrank))/real(numpartix(2),mytype)
    !   dz=(lzmax(nrank)-lzmin(nrank))/real(numpartix(3),mytype)
    !   !
    !   p=p+1
    !   !
    !   particle(p)%x=0.0001d0
    !   !
    !   particle(p)%y=dy*real(j,mytype)+lymin(nrank)-0.5d0*dy
    !   particle(p)%z=dz*real(k,mytype)+lzmin(nrank)-0.5d0*dz
    !   !
    !   particle(p)%rankinn=nrank
    !   !
    !   particle(p)%out=.false.
    !   particle(p)%swap=.false.
    !   !
    ! enddo
    ! enddo
    ! enddo
    !
    call partical_domain_check
    !
    call partical_swap
    !
    call write_particle()
    !
    part_time=0.d0
    part_comm_time=0.d0
    part_vel_time=0.d0
    part_dmck_time=0.d0
    a2a_time=0.d0
    ! call mpistop
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
  subroutine partile_inject(numadd)
    !
    use param,     only : xlx,yly,zlz
    !
    integer,intent(in) :: numadd
    !
    ! local data
    real(mytype),allocatable,dimension(:) :: xpa_add,ypa_add,zpa_add
    real(mytype),allocatable,dimension(:) :: ux_pa_add,uy_pa_add,uz_pa_add
    integer :: i,j,k,p
    !
    allocate(xpa_add(numparticle+numadd),ypa_add(numparticle+numadd),  &
             zpa_add(numparticle+numadd),ux_pa_add(numparticle+numadd),&
             uy_pa_add(numparticle+numadd),uz_pa_add(numparticle+numadd))
    !
    xpa_add(1:numparticle)=xpa(1:numparticle)
    ypa_add(1:numparticle)=ypa(1:numparticle)
    zpa_add(1:numparticle)=zpa(1:numparticle)
    !
    ux_pa_add(1:numparticle)=ux_pa(1:numparticle)
    uy_pa_add(1:numparticle)=uy_pa(1:numparticle)
    uz_pa_add(1:numparticle)=uz_pa(1:numparticle)
    !
    p=numparticle
    do k=1,numpartix(3)
    do j=1,numpartix(2)
    do i=1,numpartix(1)
      !
      p=p+1
      !
      xpa_add(p)=(partirange(2)-partirange(1))/real(numpartix(1),mytype)*  &
                 real(i,mytype)+partirange(1)
      ypa_add(p)=(partirange(4)-partirange(3))/real(numpartix(2),mytype)*  &
                  real(j,mytype)+partirange(3)
      zpa_add(p)=(partirange(6)-partirange(5))/real(numpartix(3),mytype)*  &
                 real(k,mytype)+partirange(5)
    enddo
    enddo
    enddo
    !
    deallocate(xpa,ypa,zpa,ux_pa,uy_pa,uz_pa)
    !
    call move_alloc(xpa_add, xpa)
    call move_alloc(ypa_add, ypa)
    call move_alloc(zpa_add, zpa)
    call move_alloc(ux_pa_add, ux_pa)
    call move_alloc(uy_pa_add, uy_pa)
    call move_alloc(uz_pa_add, uz_pa)
    !
    numparticle=numparticle+numadd
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
      lymax(nrank)=real(xend(2)-2,mytype)*dy
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
  !| This subroutine is to search particles.                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-11-2021  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine partivelo(ux1,uy1,uz1)
    !
    use MPI
    use param,     only : dx,dy,dz,istret,nclx,ncly,nclz,xlx,yly,zlz
    use variables, only : yp,ny,nz
    use decomp_2d, only : xsize,xstart,xend,update_halo
    use actuator_line_model_utils
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    !
    ! local data
    integer :: jpart,npart,i,j,k,psize
    type(partype),pointer :: pa
    real(mytype) :: x1,y1,z1,x2,y2,z2
    real(mytype),allocatable,dimension(:,:,:) :: ux1_halo,uy1_halo,    &
                                                 uz1_halo,ux1_hal2,    &
                                                 uy1_hal2,uz1_hal2
    !
    real(mytype),allocatable,save :: xx(:),yy(:),zz(:)
    logical,save :: firstcal=.true.
    !
    real(8) :: timebeg
    !
    timebeg=ptime()
    !
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
      do k=0,xsize(3)+1
        zz(k)=real((k+xstart(3)-2),mytype)*dz
      enddo
      !
      ! write(rankname,'(i4.4)')nrank
      ! open(18,file='testout/tecfield'//rankname//'.dat')
      ! write(18,'(A)')' VARIABLES = "x" "y" "z" "u_halo"'
      ! write(18,'(A)')'ZONE T="ZONE 001"'
      ! write(18,'(3(A,I0),A)')'I=',xsize(1),', J=',xsize(2),', K=',xsize(3),', ZONETYPE=Ordered'
      ! write(18,'(A)')'DATAPACKING=POINT'
      ! do k=1,xsize(3)
      ! do j=1,xsize(2)
      ! do i=1,xsize(1)
      !   write(18,'(4(1X,E20.13E2))')xx(i),yy(j),zz(k),ux1(i,j,k)
      ! enddo
      ! enddo
      ! enddo
    ! print*,' << testout/tecfield',rankname,'.dat'
      !
      ! print*,nrank,'|x',xmin,'~',xmax
      ! print*,nrank,'|y',ymin,'~',ymax
      ! print*,nrank,'|z',zmin,'~',real((xsize(3)+xstart(3)-1),mytype)*dz
      ! print*,nrank,'|x',xstart(1),'~',xend(1)
      ! print*,nrank,'|y',xstart(2),'~',xend(2)
      ! print*,nrank,'|z',xstart(3),'~',xend(3)
      ! print*,nrank,'|size',xsize(:),'~',size(ux1,1),size(ux1,2),size(ux1,3)
      firstcal=.false.
      !
    endif
    !
    ! use the periodic condtion, to put the particles into the domain
    ! do p=1,numparticle
    !   if(nclx .and. xpa(p)>=xlx) xpa(p)=xpa(p)-xlx
    !   if(ncly .and. ypa(p)>=yly) ypa(p)=ypa(p)-yly
    !   if(nclz .and. zpa(p)>=zlz) zpa(p)=zpa(p)-zlz
    ! enddo
    !
    call update_halo(ux1,ux1_halo,1,opt_global=.true.)
    call update_halo(uy1,uy1_halo,1,opt_global=.true.)
    call update_halo(uz1,uz1_halo,1,opt_global=.true.)
    !
    allocate(ux1_hal2(1:xsize(1),0:xsize(2)+1,0:xsize(3)+1))
    allocate(uy1_hal2(1:xsize(1),0:xsize(2)+1,0:xsize(3)+1))
    allocate(uz1_hal2(1:xsize(1),0:xsize(2)+1,0:xsize(3)+1))
    !
    ux1_hal2=ux1_halo
    uy1_hal2=uy1_halo
    uz1_hal2=uz1_halo
    !
    ! allocate(ux1_halo2(1:size(ux1_halo,1),1:size(ux1_halo,2),1:size(ux1_halo,3)))
    ! ux1_halo2=ux1_halo
    ! print*, nrank, shape(ux1), shape(ux1_halo)
    ! do k=1,xsize(3)
    ! do j=1,xsize(2)
    !   ! write(*,'(A,I0,2(I5),2(1X,E20.13E2))')'rank',nrank,j,k,ux1_halo(1,j+xstart(2)-2,k+xstart(3)-1),ux1_halo2(1,j,k)
    !   write(*,'(A,I0,2(I5),2(1X,E20.13E2))')'rank',nrank,j,k,ux1(1,j,k), &
    !                                                     ux1_halo(1,j,k)
    ! enddo
    ! enddo
    ! j=2
    ! do k=1,xsize(3)
    !   write(*,'(A,I0,2(I5),2(1X,E20.13E2))')'rank',nrank,j,k,ux1(1,j,k), &
    !                                                     ux1_hal2(1,j,k)
    ! enddo
    !
    ! print*,nrank,'|',ux1(1,:,1)
    ! print*,nrank,'|',ux1_halo(1,:,1)
    ! endif
    !
    ! write(rankname,'(i4.4)')nrank
    ! open(18,file='testout/techalo'//rankname//'.dat')
    ! write(18,'(A)')' VARIABLES = "y" "z" "u_halo"'
    ! write(18,'(A)')'ZONE T="ZONE 001"'
    ! write(18,'(2(A,I0),A)')'I=1, J=',xsize(2)+1,', K=',xsize(3)+1,', ZONETYPE=Ordered'
    ! write(18,'(A)')'DATAPACKING=POINT'
    ! do k=0,xsize(3)
    ! do j=1,xsize(2)+1
    !   write(18,'(4(1X,E20.13E2))')yy(j),zz(k),ux1_hal2(1,j,k)
    ! enddo
    ! enddo
    ! print*,' << testout/tecfield',rankname,'.dat'
    ! !
    !
    psize=msize(particle)
    !
    npart=0
    !
    do jpart=1,psize
      !
      ! print*,xpa(:),ypa(:),zpa(:)
      ! print*,nrank,'|x',(xpa(p)>=xmin .and. xpa(p)<xmax)
      ! print*,nrank,'|y',(ypa(p)>=ymin .and. ypa(p)<ymax)
      ! print*,nrank,'|z',(zpa(p)>=zmin .and. zpa(p)<zmax)
      !
      pa=>particle(jpart)
      !
      if(pa%new) cycle
      !
      loopk: do k=1,xsize(3)
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
            if( pa%x>=x1 .and. pa%x<x2 .and. &
                pa%y>=y1 .and. pa%y<y2 .and. &
                pa%z>=z1 .and. pa%z<z2 ) then
              !
              ! locate the particle, do the interpolation
              ! print*,x1,x2,y1,y2,z1,z2
              pa%u=trilinear_interpolation( x1,y1,z1,            &
                                            x2,y2,z2,            &
                                            pa%x,pa%y,pa%z,      &
                                            ux1_hal2(i,j,k),     &
                                            ux1_hal2(i+1,j,k),   &
                                            ux1_hal2(i,j,k+1),   &
                                            ux1_hal2(i+1,j,k+1), &
                                            ux1_hal2(i,j+1,k),   &
                                            ux1_hal2(i+1,j+1,k), &
                                            ux1_hal2(i,j+1,k+1), &
                                            ux1_hal2(i+1,j+1,k+1))
              pa%v=trilinear_interpolation( x1,y1,z1,            &
                                            x2,y2,z2,            &
                                            pa%x,pa%y,pa%z,      &
                                            uy1_hal2(i,j,k),     &
                                            uy1_hal2(i+1,j,k),   &
                                            uy1_hal2(i,j,k+1),   &
                                            uy1_hal2(i+1,j,k+1), &
                                            uy1_hal2(i,j+1,k),   &
                                            uy1_hal2(i+1,j+1,k), &
                                            uy1_hal2(i,j+1,k+1), &
                                            uy1_hal2(i+1,j+1,k+1)) 
              pa%w=trilinear_interpolation( x1,y1,z1,            &
                                            x2,y2,z2,            &
                                            pa%x,pa%y,pa%z,      &
                                            uz1_hal2(i,j,k),     &
                                            uz1_hal2(i+1,j,k),   &
                                            uz1_hal2(i,j,k+1),   &
                                            uz1_hal2(i+1,j,k+1), &
                                            uz1_hal2(i,j+1,k),   &
                                            uz1_hal2(i+1,j+1,k), &
                                            uz1_hal2(i,j+1,k+1), &
                                            uz1_hal2(i+1,j+1,k+1)) 
              !
              ! print*,p,ypa(p),':',ux_pa(p),uy_pa(p),uz_pa(p)
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
    ! ux_pa=psum(ux_pa)
    ! uy_pa=psum(uy_pa)
    ! uz_pa=psum(uz_pa)
    !
    ! print*,1,ypa(1),':',ux_pa(1),uy_pa(1),uz_pa(1)
    !
    ! call mpistop()
    !
  end subroutine partivelo
  !+-------------------------------------------------------------------+
  ! The end of the subroutine partivelo                                |
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
  subroutine intt_particel
    !
    use param, only : gdt,itr
    !
    ! local data 
    integer :: p,psize,jpart,npart
    type(partype),pointer :: pa
    real(8) :: timebeg
    !
    timebeg=ptime()
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
      pa%x=gdt(itr)*pa%u+pa%x
      pa%y=gdt(itr)*pa%v+pa%y
      pa%z=gdt(itr)*pa%w+pa%z
      !
      npart=npart+1
      !
      if(npart==numparticle) exit
      !
    enddo
    !
    call partical_domain_check
    !
    call partical_swap
    !
    part_time=part_time+ptime()-timebeg
    !
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
  !+-------------------------------------------------------------------+
  !| This subroutine is to check if the particle is out of domain      |
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-06-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine partical_domain_check
    !
    use param,     only : nclx,ncly,nclz,xlx,yly,zlz
    !
    ! local data 
    integer :: jpart,npart,psize,jrank,npcanc
    type(partype),pointer :: pa
    real(8) :: timebeg
    !
    timebeg=ptime()
    !
    psize=msize(particle)
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
      if(nclx .and. (pa%x>xlx .or. pa%x<0)) then
        call pa%reset()
        npcanc=npcanc+1
        cycle
      endif
      !
      if(ncly .and. (pa%y>yly .or. pa%y<0)) then
        call pa%reset()
        npcanc=npcanc+1
        cycle
      endif
      !
      if(nclz .and. (pa%z>zlz .or. pa%z<0)) then
        call pa%reset()
        npcanc=npcanc+1
        cycle
      endif
      !
      if( pa%x>=lxmin(nrank) .and. pa%x<lxmax(nrank) .and. &
          pa%y>=lymin(nrank) .and. pa%y<lymax(nrank) .and. &
          pa%z>=lzmin(nrank) .and. pa%z<lzmax(nrank) ) then
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
          if( pa%x>=lxmin(jrank) .and. pa%x<lxmax(jrank) .and. &
              pa%y>=lymin(jrank) .and. pa%y<lymax(jrank) .and. &
              pa%z>=lzmin(jrank) .and. pa%z<lzmax(jrank) ) then
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
    part_dmck_time=part_dmck_time+ptime()-timebeg
    ! print*,nrank,'|',numparticle
    ! do p=1,psize
    !   !
    !   pa=>particle(p)
    !   !
    !   if(pa%swap) then
    !     !
    !     write(*,'(3(A,1X,I0))')' ** particle',p,' moves from rank', &
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
    real(8) :: timebeg
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
    ! synchronize recv table according to send table
    nrecv=ptabupd(nsend)
    !
    ! print*,' ** rank',nrank,' total number to send',nsend_total
    !
    ! if(nrank==0) then
    !   write(*,*)' ** rank',nrank,' to send ',nsend
    !   !
    !   do jrank=0,nproc-1
    !     do jpart=1,nsend(jrank)
    !       write(*,*)' ** jrank',jrank,'|',jpart,pr(jrank,jpart)
    !     enddo
    !   enddo
    !   !
    ! endif
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
          ! write(*,*)' **',n,' particle ',p,'-> rank',jrank
          ! write(*,*)' **',n,' particle ',particle(p)%x
          !
        enddo
        !
      enddo
      !
    endif 
    !
    ! write(*,*)' ** rank',nrank,' n ',n,numparticle
    !
    ! write(*,*)' ** rank',nrank,' to recv ',nrecvtable
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! swap particle among ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call pa2a(pa2send,pa2recv,nsend,nrecv)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of swap particle among ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! now add the received particle in to the array, dynamically
    if(numparticle+msize(pa2recv)>psize) then
      !
      ! expand the particle array
      newsize=max(numparticle+msize(pa2recv),numparticle+100)
      !
      call mextend(particle,newsize)
      !
    endif
    !
    n=0
    do jpart=1,msize(particle)
      !
      pa=>particle(jpart)
      !
      ! the particle is free for re-assigning
      if(pa%new) then
        !
        if(n>=msize(pa2recv)) exit
        !
        n=n+1
        !
        pa=pa2recv(n)
        pa%new=.false.
        !
      endif
      !
    enddo
    !
    numparticle=numparticle+n
    !
    part_comm_time=part_comm_time+ptime()-timebeg
    !
    ! print*,nrank,'|',numparticle
    ! !
    ! call mpistop
    !
  end subroutine partical_swap
  !+-------------------------------------------------------------------+
  ! The end of the subroutine partical_swap                            |
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
    integer :: p,psize
    logical,save :: firstcal=.true.
    !
    if(firstcal) then
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
        write(18,'(3(1X,E15.7E3),1X,I0)')particle(p)%x,particle(p)%y,particle(p)%z,particle(p)%rankinn
      enddo
      close(18)
      print*,' << ./data/particle',rankname,'.dat'
    endif
    !
  end subroutine write_particle
  !+-------------------------------------------------------------------+
  ! The end of the subroutine write_particle                           |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is used to finalise mpi and stop the program.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-July-2019: Created by J. Fang @ STFC Daresbury Laboratory      |
  !+-------------------------------------------------------------------+
  subroutine mpistop
    !
    use mpi
    !
    integer :: ierr
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    call mpi_finalize(ierr)
    !
    if(nrank==0) print*,' ** The job is done!'
    !
    stop
    !
  end subroutine mpistop
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mpistop.                                |
  !+-------------------------------------------------------------------+
  !!
  function psum_mytype_ary(var) result(varsum)
    !
    use mpi
    use decomp_2d, only : real_type
    !
    ! arguments
    real(mytype),intent(in) :: var(:)
    real(mytype),allocatable :: varsum(:)
    !
    ! local data
    integer :: ierr,nsize
    !
    nsize=size(var)
    !
    allocate(varsum(nsize))
    !
    call mpi_allreduce(var,varsum,nsize,real_type,mpi_sum,             &
                                                    mpi_comm_world,ierr)
    !
    return
    !
  end function psum_mytype_ary
  !
  function psum_integer(var) result(varsum)
    !
    use mpi
    use decomp_2d, only : real_type
    !
    ! arguments
    integer,intent(in) :: var
    integer :: varsum
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,varsum,1,mpi_integer,mpi_sum,           &
                                                    mpi_comm_world,ierr)
    !
    return
    !
  end function psum_integer
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
  !!
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
  !| this function is to update table based on alltoall mpi            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Jun-2022  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  function ptable_update_int_arr(vain) result(vout)
    !
    use mpi
    !
    integer,intent(in) :: vain(:)
    integer :: vout(size(vain))
    !
    ! local variables
    integer :: nvar,ierr
    !
    nvar=size(vain)
    !
    call mpi_alltoall(vain,1,mpi_integer,                   &
                      vout,1,mpi_integer,mpi_comm_world,ierr)
    !
    return
    !
  end function ptable_update_int_arr
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ptable_update_int_arr.                  |
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
    real(8),allocatable :: r8send(:,:),r8resv(:,:)
    !
    integer,save :: newtype
    !
    logical,save :: firstcal=.true.
    !
    real(8) :: timebeg
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
      r8send(1,jpart)=datasend(jpart)%x
      r8send(2,jpart)=datasend(jpart)%y
      r8send(3,jpart)=datasend(jpart)%z
      r8send(4,jpart)=datasend(jpart)%u
      r8send(5,jpart)=datasend(jpart)%v
      r8send(6,jpart)=datasend(jpart)%w
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
        datarecv(jc)%x=r8resv(1,jc)
        datarecv(jc)%y=r8resv(2,jc)
        datarecv(jc)%z=r8resv(3,jc)
        datarecv(jc)%u=r8resv(4,jc)
        datarecv(jc)%v=r8resv(5,jc)
        datarecv(jc)%w=r8resv(6,jc)
        !
        ! print*,nrank,'|',datarecv(jc)%x,datarecv(jc)%y,datarecv(jc)%z
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
  !+-------------------------------------------------------------------+
  !| The wraper of MPI_Wtime                                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-November-2019: Created by J. Fang @ STFC Daresbury Laboratory  |
  !+-------------------------------------------------------------------+
  real(8) function ptime()
    !
    use mpi
    !
    ptime=MPI_Wtime()
    !
    return
    !
  end function ptime
  !+-------------------------------------------------------------------+
  !| The end of the function ptime.                                    |
  !+-------------------------------------------------------------------+
  !
end module partack
!+---------------------------------------------------------------------+
! The end of the module partack                                        |
!+---------------------------------------------------------------------+