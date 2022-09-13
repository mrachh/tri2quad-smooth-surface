      subroutine trimesh_to_quadmesh(fname_base,norder,iref,fnameout)
!
!   Given a go3 file of the form fname_base_o(norder)_r(iref), convert 
!   it to a quad mesh assuming that the at refinement level0, the 
!   triangles 2*i-1,2*i form quad i
!
!

      implicit real *8 (a-h,o-z)
      character(len=*) fname_base,fnameout
      integer iref,norder
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),wts(:)
      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      real *8, allocatable :: uvs_quad(:,:),src_quad(:,:,:)
      real *8, allocatable :: srccoefs_quad(:,:,:)
      real *8, allocatable :: u(:,:),v(:,:),w(:)
      integer, allocatable :: itri0(:,:),itri1(:,:)
      character *200 fname1,fname0

      call prini(6,13)
!
!  open geometry info
!
!
      write(fname1,'(a,a,i2.2,a,i1,a)') trim(fname_base),'_o', &
        norder,'_r0',iref,'.go3'

      call open_gov3_geometry_mem(fname1,npatches,npts)

      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,&
         iptype,npts,srcvals,srccoefs,wts)
      nquad = npatches/2
!
! figure out itri
!
!
      write(fname0,'(a,a,i2.2,a,i1,a)') trim(fname_base),'_o', &
        norder,'_r0',0,'.go3'
      call open_gov3_geometry_mem(fname0,npatches0,npts0)
      nquad0 = npatches0/2

      allocate(itri0(2,nquad),itri1(2,nquad))

      do i=1,nquad0
        itri0(1,i) = 2*i-1
        itri0(2,i) = 2*i 
      enddo

      nq = nquad0
      nq4 = 4*nq
      do i=1,iref
        call get_quad2tri_indices_ref(nq,itri0,nq4,itri1)
        itri0(1:2,1:nq4) = itri1(1:2,1:nq4)
        nq = nq4
        nq4 = 4*nq
      enddo

!
! Now convert to quads
!
!
      nqpols = (norder+1)*(norder+1)
      ntpols = (norder+1)*(norder+2)/2
      allocate(src_quad(3,nqpols,nquad),srccoefs_quad(3,nqpols,nquad))
      allocate(uvs_quad(2,nqpols),w(nqpols),u(nqpols,nqpols))
      allocate(v(nqpols,nqpols))

      itype = 2
      ipoly = 1
      call polytens_exps_2d(ipoly,itype,norder+1,uvs_quad,'f',u,nqpols, &
        v,nqpols,w)
      call prin2('uvs_quad=*',uvs_quad,24)
      
      open(unit=33,file=trim(fnameout))
      alpha = 1.0d0
      beta = 0.0d0
      do i=1,nquad
        i1 = itri0(1,i)
        i2 = itri0(2,i)
        call convert_tri_to_quad(norder,ntpols,nqpols, &
          uvs_quad,srccoefs(1,ixyzs(i1)), &
          srccoefs(1,ixyzs(i2)),src_quad(1,1,i))
        do j=1,nqpols
          write(33,'(3(2x,e22.16))') src_quad(1,j,i), &
            src_quad(2,j,i),src_quad(3,j,i)
        enddo

        call dgemm_guru('n','t',3,nqpols,nqpols,alpha,src_quad(1,1,i),3,&
          u,nqpols,beta,srccoefs_quad(1,1,i))
        if(i.le.3) call prin2('srccoefs_quad=*',srccoefs_quad(1,1,i),3*nqpols)
      enddo


      close(33)

      return
      end
!
!
!
!
!




      subroutine convert_tri_to_quad(norder,ntpols,nqpols, &
        uvs_quad,srccoefs_tri1,srccoefs_tri2,srcinfo_quad)
!
!  This subroutine converts Koornwinder polynomial expansions
!  of two triangles which form a quad into a a chebyshev expansion for
!  a corresponding quad. The quad is assumed to be discretized on
!  (-1,1)^2
!
!
!  Input arguments:
!    - norder: integer
!        order of discretization
!    - ntpols: integer
!        number of polynomials on the triangle
!    - nqpols: integer
!        number of polynomials on the quad
!    - uvs_quad: real *8 (2,nqpols)
!        u,v coordinates of (0,1)^2 quad
!    - srccoefs_tri1: real *8 (9,ntpols)
!        Koornwinder expansion coefs of triangle 1
!    - srccoefs_tri2: real *8 (9,ntpols)
!        Koornwinder expansion coefs of triangle 2
!
!  Output arguments:
!    - srcinfo_quad: real *8 (3,nqpols)
!        xyz values of quad
!
!  Todo: 
!    Figure out how to compute derivatives correctly
!
      implicit real *8 (a-h,o-z)
      integer norder,ntpols,nqpols
      real *8 srccoefs_tri1(9,ntpols),srccoefs_tri2(9,ntpols)
      real *8 uvs_quad(2,nqpols)
      real *8 srcinfo_quad(3,nqpols),pols(ntpols),uvs_targ(2)
      

      do i=1,nqpols
        u = uvs_quad(1,nqpols)
        v = uvs_quad(2,nqpols)
        srcinfo_quad(1:3,i) = 0

        if(u+v.le.1) then 
          uvs_targ(1) = u
          uvs_targ(2) = v
          call koorn_pols(uvs_targ,norder,ntpols,pols)
          do j=1,ntpols
            srcinfo_quad(1:3,i) = srcinfo_quad(1:3,i) + &
              srccoefs_tri1(1:3,j)*pols(j)
          enddo
        else
          uvs_targ(1) = 1.0d0-u
          uvs_targ(2) = 1.0d0-v
          call koorn_pols(uvs_targ,norder,ntpols,pols)
          do j=1,ntpols
            srcinfo_quad(1:3,i) = srcinfo_quad(1:3,i) + &
              srccoefs_tri2(1:3,j)*pols(j)
          enddo
        endif
      enddo


      return
      end


      subroutine get_quad2tri_indices_ref(nq,itri,nq4,itri4)
!
!  If a quad \tilde{i} consists of triangles i,j, stored in itri, then
!  which children of i,j form the children of quad \tilde{i} are
!  stored in itriref
!   
!  Input arguments:
!    - nq: integer
!        number of quads
!    - itri: integer(2,nq)
!        which triangles form quads i
!    - nq4: integer
!        number of children quad
!
!  Output arguments:
!    - itri4: integer (2,nq4)
!        which triangles for the children of quad \tilde{i}
!
!  Notes: the children of triangle i in the refined list are:
!    4i-3,4i-2,4i-1, and 4i. Similarly, the children of quad i
!    are 4i-3,4i-2,4i-1,4i
!
!
      implicit real *8 (a-h,o-z)
      integer nq,nq4,itri(2,nq),itri4(2,nq4)

      do i=1,nq
        itri4(1,4*i-3) = 4*itri(1,i)-3
        itri4(2,4*i-3) = 4*itri(1,i)-2
        
        itri4(1,4*i-2) = 4*itri(1,i)-1
        itri4(2,4*i-2) = 4*itri(2,i)
        
        itri4(1,4*i-1) = 4*itri(1,i)
        itri4(2,4*i-1) = 4*itri(2,i)-1

        itri4(1,4*i) = 4*itri(2,i)-2
        itri4(2,4*i) = 4*itri(2,i)-3
      enddo

      return
      end
