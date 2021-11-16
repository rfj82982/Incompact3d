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
  interface psum
    module procedure psum_mytype_ary
  end interface
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
    use decomp_2d, only : xsize,xstart,xend,nrank,update_halo
    use actuator_line_model_utils
    !
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    !
    ! local data
    integer :: p,nyr,nzr,i,j,k
    real(mytype) :: x1,y1,z1,x2,y2,z2
    real(mytype),allocatable,dimension(:,:,:) :: ux1_halo,uy1_halo,    &
                                                 uz1_halo,ux1_hal2,    &
                                                 uy1_hal2,uz1_hal2
    real(mytype),save :: xmin,xmax,ymin,ymax,zmin,zmax
    real(mytype),allocatable,save :: xx(:),yy(:),zz(:)
    logical,save :: firstcal=.true.
    !
    character(len=4) :: rankname
    !
    if(firstcal) then
      !
      allocate(xx(0:xsize(1)+1),yy(0:xsize(2)+1),zz(0:xsize(3)+1))
      !
      xmin=real(0,mytype)*dx
      xmax=real(xend(1),mytype)*dx
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
      if (istret==0) then
        ymin=real(xstart(2)-1,mytype)*dy
        ymax=real(xend(2)-2,mytype)*dy
      else
        ymin=yp(xstart(2))
        nyr=min(ny,xend(2)+1)
        ymax=yp(nyr)
      endif
      !
      zmin=real((xstart(3)-1),mytype)*dz
      nzr=xend(3)
      zmax=real(nzr,mytype)*dz
      !
      write(rankname,'(i4.4)')nrank
      open(18,file='testout/tecfield'//rankname//'.dat')
      write(18,'(A)')' VARIABLES = "x" "y" "z" "u_halo"'
      write(18,'(A)')'ZONE T="ZONE 001"'
      write(18,'(3(A,I0),A)')'I=',xsize(1),', J=',xsize(2),', K=',xsize(3),', ZONETYPE=Ordered'
      write(18,'(A)')'DATAPACKING=POINT'
      do k=1,xsize(3)
      do j=1,xsize(2)
      do i=1,xsize(1)
        write(18,'(4(1X,E20.13E2))')xx(i),yy(j),zz(k),ux1(i,j,k)
      enddo
      enddo
      enddo
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
    do p=1,numparticle
      if(nclx .and. xpa(p)>=xlx) xpa(p)=xpa(p)-xlx
      if(ncly .and. ypa(p)>=yly) ypa(p)=ypa(p)-yly
      if(nclz .and. zpa(p)>=zlz) zpa(p)=zpa(p)-zlz
    enddo
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
    do p=1,numparticle
      !
      ! print*,xpa(:),ypa(:),zpa(:)
      ! print*,nrank,'|x',(xpa(p)>=xmin .and. xpa(p)<xmax)
      ! print*,nrank,'|y',(ypa(p)>=ymin .and. ypa(p)<ymax)
      ! print*,nrank,'|z',(zpa(p)>=zmin .and. zpa(p)<zmax)
      !
      if( xpa(p)>=xmin .and. xpa(p)<xmax .and. &
          ypa(p)>=ymin .and. ypa(p)<ymax .and. &
          zpa(p)>=zmin .and. zpa(p)<zmax ) then
        !
        loopk: do k=1,xsize(3)
          z1=zz(k)
          z2=zz(k+1)
          !
          ! print*,k,z1,z2
          !
          do j=1,xsize(2)
            y1=yy(j)
            y2=yy(j+1)
            do i=1,xsize(1)
              x1=xx(i)
              x2=xx(i+1)
              !
              if( xpa(p)>=x1 .and. xpa(p)<x2 .and. &
                  ypa(p)>=y1 .and. ypa(p)<y2 .and. &
                  zpa(p)>=z1 .and. zpa(p)<z2 ) then
                !
                ! locate the particle, do the interpolation
                ! print*,x1,x2,y1,y2,z1,z2
                ux_pa(p)=trilinear_interpolation( x1,y1,z1, &
                                                  x2,y2,z2, &
                                                  xpa(p),ypa(p),zpa(p),&
                                                  ux1_hal2(i,j,k),     &
                                                  ux1_hal2(i+1,j,k),   &
                                                  ux1_hal2(i,j,k+1),   &
                                                  ux1_hal2(i+1,j,k+1), &
                                                  ux1_hal2(i,j+1,k),   &
                                                  ux1_hal2(i+1,j+1,k), &
                                                  ux1_hal2(i,j+1,k+1), &
                                                  ux1_hal2(i+1,j+1,k+1))
                uy_pa(p)=trilinear_interpolation( x1,y1,z1, &
                                                  x2,y2,z2, &
                                                  xpa(p),ypa(p),zpa(p),&
                                                  uy1_hal2(i,j,k),     &
                                                  uy1_hal2(i+1,j,k),   &
                                                  uy1_hal2(i,j,k+1),   &
                                                  uy1_hal2(i+1,j,k+1), &
                                                  uy1_hal2(i,j+1,k),   &
                                                  uy1_hal2(i+1,j+1,k), &
                                                  uy1_hal2(i,j+1,k+1), &
                                                  uy1_hal2(i+1,j+1,k+1)) 
                uz_pa(p)=trilinear_interpolation( x1,y1,z1, &
                                                  x2,y2,z2, &
                                                  xpa(p),ypa(p),zpa(p),&
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
      else
        !
        ux_pa(p)=0.d0
        uy_pa(p)=0.d0
        uz_pa(p)=0.d0
        !
      endif 
      !
    enddo
    !
    ux_pa=psum(ux_pa)
    uy_pa=psum(uy_pa)
    uz_pa=psum(uz_pa)
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
    use decomp_2d, only : nrank
    use param, only : gdt,itr
    !
    ! local data 
    integer :: p
    !
    do p=1,numparticle
      xpa(p)=gdt(itr)*ux_pa(p)+xpa(p)
      ypa(p)=gdt(itr)*uy_pa(p)+ypa(p)
      zpa(p)=gdt(itr)*uz_pa(p)+zpa(p)
    enddo
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
    use decomp_2d, only : nrank
    !
    ! local data
    integer :: p
    logical,save :: firstcal=.true.
    !
    if(nrank==0) then
      !
      if(firstcal) then
        open(18,file='./data/particle.dat')
        write(18,'(A)')'VARIABLES = "x" "y" "z"'
        close(18)
        print*,' create particle.dat'
        !
        firstcal=.false.
      endif
      !
      open(18,file='./data/particle.dat',position="append")
      write(18,'(A)')'ZONE T="ZONE 001"'
      write(18,'(A,I0,A)')'I=',numparticle,', J=1, K=1, ZONETYPE=Ordered'
      write(18,'(A)')'DATAPACKING=POINT'
      do p=1,numparticle
        write(18,*)xpa(p),ypa(p),zpa(p)
      enddo
      close(18)
      print*,' << ./data/particle.dat'
      !
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
    use decomp_2d, only : nrank
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
end module particle
!+---------------------------------------------------------------------+
! The end of the module particle                                       |
!+---------------------------------------------------------------------+