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
  use hdf5
  use h5lt
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
    module procedure updatable_int
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
  Interface h5write
    !
    module procedure h5wa_r8
    module procedure h5w_real8
    module procedure h5w_int4
    !
  end Interface h5write
  !
  interface pgather
    module procedure pgather_int
  end interface
  ! particles
  type partype
    !
    real(mytype) :: x(3),v(3)
    real(mytype),allocatable,dimension(:,:) :: dx,dv
    integer :: id,rankinn,rank2go
    logical :: swap,new
    !+------------------+------------------------------------------+
    !|                x | spatial coordinates of particle          |
    !|                v | velocity of particle                     |
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
  real(8) :: part_time,part_comm_time,part_vel_time,part_dmck_time,a2a_time, &
             count_time,data_pack_time,data_unpack_time,mpi_comm_time,       &
             table_share_time
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
  integer(hid_t) :: h5file_id
  integer :: mpi_comm_particle,mpi_rank_part,mpi_size_part
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
  end subroutine init_one_particle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine initmesg.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to reset a particle.                       |
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
    real(8) :: dx,dy,dz,x,y,z
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
      if( x>=lxmin(nrank) .and. x<lxmax(nrank) .and. &
          y>=lymin(nrank) .and. y<lymax(nrank) .and. &
          z>=lzmin(nrank) .and. z<lzmax(nrank) ) then
        !
        p=p+1
        !
        call particle_new(p)%init()
        !
        particle_new(p)%x(1)=x
        particle_new(p)%x(2)=y
        particle_new(p)%x(3)=z
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
    use param,     only : xlx,yly,zlz
    use var,       only : itime
    !
    ! local data
    integer :: i,j,k,p
    real(mytype) :: dx,dy,dz
    !
    call particle_gen(particle,numparticle)
    !
    ! call partical_domain_check
    ! !
    ! call partical_swap
    !
    ! call write_particle()
    call h5write_particle(0)
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
    call partical_domain_check
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
  subroutine particle_velo(ux1,uy1,uz1)
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
            if( pa%x(1)>=x1 .and. pa%x(1)<x2 .and. &
                pa%x(2)>=y1 .and. pa%x(2)<y2 .and. &
                pa%x(3)>=z1 .and. pa%x(3)<z2 ) then
              !
              ! locate the particle, do the interpolation
              ! print*,x1,x2,y1,y2,z1,z2
              pa%v(1)=trilinear_interpolation( x1,y1,z1,            &
                                            x2,y2,z2,            &
                                            pa%x(1),pa%x(2),pa%x(3),      &
                                            ux1_hal2(i,j,k),     &
                                            ux1_hal2(i+1,j,k),   &
                                            ux1_hal2(i,j,k+1),   &
                                            ux1_hal2(i+1,j,k+1), &
                                            ux1_hal2(i,j+1,k),   &
                                            ux1_hal2(i+1,j+1,k), &
                                            ux1_hal2(i,j+1,k+1), &
                                            ux1_hal2(i+1,j+1,k+1))
              pa%v(2)=trilinear_interpolation( x1,y1,z1,            &
                                            x2,y2,z2,            &
                                            pa%x(1),pa%x(2),pa%x(3),      &
                                            uy1_hal2(i,j,k),     &
                                            uy1_hal2(i+1,j,k),   &
                                            uy1_hal2(i,j,k+1),   &
                                            uy1_hal2(i+1,j,k+1), &
                                            uy1_hal2(i,j+1,k),   &
                                            uy1_hal2(i+1,j+1,k), &
                                            uy1_hal2(i,j+1,k+1), &
                                            uy1_hal2(i+1,j+1,k+1)) 
              pa%v(3)=trilinear_interpolation( x1,y1,z1,            &
                                            x2,y2,z2,            &
                                            pa%x(1),pa%x(2),pa%x(3),      &
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
  end subroutine particle_velo
  !+-------------------------------------------------------------------+
  ! The end of the subroutine particle_velo                            |
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
    use variables 
    use param
    !
    ! local data 
    integer :: p,psize,jpart,npart,total_num_part
    integer,save :: old_num_part=0
    type(partype),pointer :: pa
    real(8) :: timebeg
    !
    real(8),allocatable :: xcor(:,:),dxco(:,:,:)
    !
    timebeg=ptime()
    !
    allocate(xcor(3,numparticle),dxco(3,numparticle,ntime))
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
      dxco(:,npart,1)=pa%v
      dxco(:,npart,2:ntime)=pa%dx(:,2:ntime)
      !
      if(npart==numparticle) exit
      !
    enddo
    !
    if (itimescheme.eq.1) then
       !>>> Euler
       xcor=gdt(itr)*dxco(:,:,1)+xcor
       !
    elseif(itimescheme.eq.2) then
       !>>> Adam-Bashforth second order (AB2)
       !
       if(itime.eq.1 .and. irestart.eq.0) then
         ! Do first time step with Euler
         xcor=gdt(itr)*dxco(:,:,1)+xcor
       else
         xcor=adt(itr)*dxco(:,:,1)+bdt(itr)*dxco(:,:,2)+xcor
       endif
       dxco(:,:,2)=dxco(:,:,1)
       !
    elseif(itimescheme.eq.3) then
       !>>> Adams-Bashforth third order (AB3)
       !
       ! Do first time step with Euler
       if(itime.eq.1.and.irestart.eq.0) then
          xcor=dt*dxco(:,:,1)+xcor
       elseif(itime.eq.2.and.irestart.eq.0) then
          ! Do second time step with AB2
          xcor=onepfive*dt*dxco(:,:,1)-half*dt*dxco(:,:,2)+xcor
          dxco(:,:,3)=dxco(:,:,2)
       else
          ! Finally using AB3
          xcor=adt(itr)*dxco(:,:,1)+bdt(itr)*dxco(:,:,2)+cdt(itr)*dxco(:,:,3)+xcor
          dxco(:,:,3)=dxco(:,:,2)
       endif
       dxco(:,:,2)=dxco(:,:,1)
       !
    elseif(itimescheme.eq.5) then
       !>>> Runge-Kutta (low storage) RK3
       if(itr.eq.1) then
          xcor=gdt(itr)*dxco(:,:,1)+xcor
       else
          xcor=adt(itr)*dxco(:,:,1)+bdt(itr)*dxco(:,:,2)+xcor
       endif
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
      pa%x(:)=xcor(:,npart)
      !
      pa%dx(:,2:ntime)=dxco(:,npart,2:ntime)
      !
      if(npart==numparticle) exit
      !
    enddo
    !
    deallocate(xcor,dxco)
    !
    call partical_domain_check
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
    integer :: jpart,npart,psize,jrank,npcanc,npcanc_totl
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
    real(8) :: timebeg,tvar1,tvar11,tvar2,tvar3,tvar4
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
  !+-------------------------------------------------------------------+
  !| This subroutine is to write particles via HDF5.                   |
  !+-------------------------------------------------------------------+
  !| tecplot format for now                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-09-2022  | Created by J. Fang                                  |
  !+-------------------------------------------------------------------+
  subroutine h5write_particle(itime)
    !
    use param, only : t
    !
    ! arguments
    integer, intent(in) :: itime
    !
    ! local data
    character(len=3) :: num
    real(8),allocatable :: xpart(:),ypart(:),zpart(:)
    integer :: psize,total_num_part,p,j
    integer :: rank2coll
    character(len=32) :: file2write
    !
    if(numparticle>0) then
      !
      allocate(xpart(numparticle),ypart(numparticle),zpart(numparticle))
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
    if(rank2coll>=0) then
      !
      write(num,'(I3.3)') itime
      !
      file2write='./data/particle'//num//'.h5'
      call h5io_init(filename=trim(file2write),mode='writ',comm=mpi_comm_particle)
      !
      call h5write(varname='time',var=t)
      call h5write(varname='x',var=xpart,total_size=total_num_part,comm=mpi_comm_particle)
      call h5write(varname='y',var=ypart,total_size=total_num_part,comm=mpi_comm_particle)
      call h5write(varname='z',var=zpart,total_size=total_num_part,comm=mpi_comm_particle)
      !
      deallocate(xpart,ypart,zpart)
      !
      call h5io_end
      !
    endif
    !
    if(nrank==0) print*,' ** total number of particles is:',total_num_part
    !
  end subroutine h5write_particle
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
  function psum_integer(var,comm) result(varsum)
    !
    use mpi
    use decomp_2d, only : real_type
    !
    ! arguments
    integer,intent(in) :: var
    integer,optional,intent(in) :: comm
    integer :: varsum
    !
    ! local data
    integer :: ierr,comm_2_use
    !
    if(present(comm)) then
        comm_2_use=comm
    else
        comm_2_use=mpi_comm_world
    endif
    !
    !
    call mpi_allreduce(var,varsum,1,mpi_integer,mpi_sum,           &
                                                    comm_2_use,ierr)
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
  !
  function updatable_int(var,offset,debug,comm) result(table)
    !
    use mpi
    !
    ! arguments
    integer :: table(0:mpi_size_part-1)
    integer,intent(in) :: var
    integer,optional,intent(out) :: offset
    logical,intent(in),optional :: debug
    integer,intent(in),optional :: comm
    !
    ! local data
    integer :: comm_2_use
    integer :: ierr,i
    integer :: vta(0:mpi_size_part-1)
    logical :: ldebug
    !
    if(present(debug)) then
      ldebug=debug
    else
      ldebug=.false.
    endif
    !
    if(present(comm)) then
        comm_2_use=comm
    else
        comm_2_use=mpi_comm_world
    endif
    !
    call mpi_allgather(var,1,mpi_integer,                              &
                       vta,1,mpi_integer,comm_2_use,ierr)
    !
    table=vta
    !
    if(present(offset)) then
      !
      if(nrank==0) then
        offset=0
      else
        !
        offset=0
        do i=0,nrank-1
          offset=offset+vta(i)
        enddo
        !
      endif
      !
    endif
    !
  end function updatable_int
  !
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
  !+-------------------------------------------------------------------+
  !| The wraper of MPI_Wtime                                           |
  !+-------------------------f------------------------------------------+
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
  !+-------------------------------------------------------------------+
  !| This subroutine is used to open the h5file interface and assign   |
  !| h5file_id. For write each new file, this will be called first, but|
  !| once it is called, the file will be overwriten.                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine h5io_init(filename,mode,comm)
    !
    use mpi, only: mpi_comm_world,mpi_info_null
    !
    ! arguments
    character(len=*),intent(in) :: filename
    character(len=*),intent(in) :: mode
    integer,intent(in),optional :: comm
    ! h5file_id is returned
    !
    ! local data
    integer :: h5error,comm_2_use
    integer(hid_t) :: plist_id
    !
    if(present(comm)) then
        comm_2_use=comm
    else
        comm_2_use=mpi_comm_world
    endif
    !
    call h5open_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_init call h5open_f'
    !
    ! create access property list and set mpi i/o
    call h5pcreate_f(h5p_file_access_f,plist_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_init call h5pcreate_f'
    !
    call h5pset_fapl_mpio_f(plist_id,comm_2_use,mpi_info_null,     &
                                                                h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_init call h5pset_fapl_mpio_f'
    !
    if(mode=='writ') then
      call h5fcreate_f(filename,h5f_acc_trunc_f,h5file_id,             &
                                            h5error,access_prp=plist_id)
      if(h5error.ne.0)  stop ' !! error in h5io_init call h5fcreate_f'
    elseif(mode=='read') then
      call h5fopen_f(filename,h5f_acc_rdwr_f,h5file_id,                &
                                            h5error,access_prp=plist_id)
      if(h5error.ne.0)  stop ' !! error in h5io_init call h5fopen_f'
    else
        stop ' !! mode not defined @ h5io_init'
    endif
    !
    call h5pclose_f(plist_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_init call h5pclose_f'
    !
    if(nrank==0) print*,' ** open h5 file: ',filename
    !
  end subroutine h5io_init
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5io_init.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to close hdf5 interface after finish      |
  !| input/output a hdf5 file.                                         |
  !| the only data needed is h5file_id                                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine h5io_end
    !
    ! local data
    integer :: h5error
    !
    call h5fclose_f(h5file_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_end call h5fclose_f'
    !
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_end call h5close_f'
    !
  end subroutine h5io_end
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5io_end.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write a 1D array with hdf5 interface.  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine h5wa_r8(varname,var,total_size,comm)
    !
    use decomp_2d, only : mytype
    use mpi, only: mpi_comm_world,mpi_info_null
    ! use parallel,only: nrank,mpistop,psum,ptabupd,nrankmax
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(mytype),intent(in),allocatable :: var(:)
    integer,intent(in) :: total_size
    integer,intent(in),optional :: comm
    !
    ! local data
    integer :: jrk,comm_2_use
    integer :: dim,dima
    integer :: dim_table(0:mpi_size_part-1)
    integer(hsize_t), dimension(1) :: offset
    integer :: h5error
    !
    integer(hid_t) :: dset_id,filespace,memspace,plist_id
    integer(hsize_t) :: dimt(1),dimat(1)
    !
    if(allocated(var)) then
      dim=size(var)
    else
      dim=0
    endif
    !
    if(present(comm)) then
        comm_2_use=comm
    else
        comm_2_use=mpi_comm_world
    endif
    !
    ! dima=psum(dim,comm=comm_2_use)
    !
    dimt=(/dim/)
    dimat=(/total_size/)
    !
    dim_table=ptabupd(dim,comm=mpi_comm_particle)
    !
    offset=0
    do jrk=0,mpi_rank_part-1
      offset=offset+dim_table(jrk)
    enddo
    !
    ! print*,mpi_rank_part,nrank,'|',dim,offset
    !
    ! writing the data
    !
    call h5screate_simple_f(1,dimat,filespace,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5screate_simple_f'
    call h5dcreate_f(h5file_id,varname,h5t_native_double,filespace,    &
                                                       dset_id,h5error)

    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5dcreate_f'
    call h5screate_simple_f(1,dimt,memspace,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5screate_simple_f'
    call h5sclose_f(filespace,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5sclose_f'
    call h5dget_space_f(dset_id,filespace,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5dget_space_f'
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,      &
                                                          dimt,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5sselect_hyperslab_f'
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5pcreate_f'
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5pset_dxpl_mpio_f'
    call h5dwrite_f(dset_id,h5t_native_double,var,dimt,h5error,        &
                    file_space_id=filespace,mem_space_id=memspace,     &
                                                     xfer_prp=plist_id)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5dwrite_f'
    call h5sclose_f(filespace,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5sclose_f'
    call h5sclose_f(memspace,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5sclose_f'
    call h5dclose_f(dset_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5dclose_f'
    call h5pclose_f(plist_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5wa_r8 call h5pclose_f'
    !
    if(mpi_rank_part==0) print*,' << ',varname
    !
  end subroutine h5wa_r8
  !
  subroutine h5w_int4(varname,var)
    !
    use mpi, only: mpi_info_null
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    integer,intent(in) :: var
    !
    ! local data
    integer :: nvar(1)
    integer :: h5error
    integer(hsize_t) :: dimt(1)=(/1/)
    !
    ! writing the data
    !
    nvar=var
    call h5ltmake_dataset_f(h5file_id,varname,1,dimt,                  &
                                        h5t_native_integer,nvar,h5error)
    if(h5error.ne.0)  stop ' !! error in h5w_int4 call h5ltmake_dataset_f'
    !
    if(nrank==0) print*,' << ',varname
    !
  end subroutine h5w_int4
  !
  subroutine h5w_real8(varname,var)
    !
    use mpi, only: mpi_info_null
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(in) :: var
    !
    ! local data
    real(8) :: rvar(1)
    integer :: h5error
    integer(hsize_t) :: dimt(1)=(/1/)
    !
    rvar=var
    call h5ltmake_dataset_f(h5file_id,varname,1,dimt,                  &
                                        h5t_native_double,rvar,h5error)
    if(h5error.ne.0)  stop ' !! error in h5w_real8 call h5ltmake_dataset_f'
    !
    if(nrank==0) print*,' << ',varname
    !
  end subroutine h5w_real8
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to create a sub-communicator from nranks.    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-08-2022: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine subcomm_group(rank,communicator,newrank,newsize)
    !
    use mpi
    ! arguments
    integer,intent(in) :: rank
    integer,intent(out) :: communicator,newrank,newsize
    !
    ! local data
    integer :: group_mpi,mpi_group_world
    integer :: ierr,ncout,jrank
    integer,allocatable :: rank_use(:),ranktemp(:)
    !
    allocate(ranktemp(0:nproc-1))
    !
    call pgather(rank,ranktemp)
    !
    ncout=0
    do jrank=0,nproc-1
      !
      if(ranktemp(jrank)>=0) then
        ncout=ncout+1
      endif
      !
    enddo
    !
    allocate(rank_use(1:ncout))
    !
    ncout=0
    do jrank=0,nproc-1
      !
      if(ranktemp(jrank)>=0) then
        ncout=ncout+1
        !
        rank_use(ncout)=ranktemp(jrank)
        !
      endif
      !
    enddo
    !
    call mpi_comm_group(mpi_comm_world,mpi_group_world,ierr)
    call mpi_group_incl(mpi_group_world,size(rank_use),rank_use,group_mpi,ierr)
    call mpi_comm_create(mpi_comm_world,group_mpi,communicator,ierr)
    !
    if(any(rank_use==nrank)) then
      call mpi_comm_size(communicator,newsize,ierr)
      call mpi_comm_rank(communicator,newrank,ierr)
      if(newrank==0) print*,' ** new subcomm created, size: ',newsize
      ! print*,' ** local rank:',newrank,', gloable rank:',nrank
    else
      newrank=-1
      newsize=0
    endif
    !
  end subroutine subcomm_group
  !+-------------------------------------------------------------------+
  !| The end of the subroutine subcomm_group.                          |
  !+-------------------------------------------------------------------+
  !
  subroutine pgather_int(var,data,mode)
    !
    use mpi
    !
    ! arguments
    integer,intent(in) :: var
    integer,intent(out),allocatable :: data(:)
    character(len=*),intent(in),optional :: mode
    !
    !
    ! local data
    integer :: counts(0:nproc-1)
    integer :: ierr,jrank,ncou
    !
    call mpi_allgather(var, 1, mpi_integer, counts, 1, mpi_integer,  &
                       mpi_comm_world, ierr)
    !
    if(present(mode) .and. mode=='noneg') then
      ! only pick >=0 values
      ncou=0
      do jrank=0,nproc-1
        if(counts(jrank)>=0) then
          ncou=ncou+1
        endif
      enddo
      !
      allocate(data(ncou))
      ncou=0
      do jrank=0,nproc-1
        if(counts(jrank)>=0) then
          ncou=ncou+1
          data(ncou)=counts(jrank)
        endif
      enddo
      !
    else
      allocate(data(0:nproc-1))
      data=counts
    endif
    !
  end subroutine pgather_int
  !
end module partack
!+---------------------------------------------------------------------+
! The end of the module partack                                        |
!+---------------------------------------------------------------------+