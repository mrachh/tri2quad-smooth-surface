      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:)
      complex *16, allocatable :: sigma2(:),rhs2(:)

      complex ( kind = 8 ), allocatable :: a_vect(:,:),RHS_vect(:)
      real *8, allocatable :: dplot(:,:)
      real *8, allocatable :: ru(:,:),rv(:,:)
      complex *16 vf(3)

      integer nnn,mmm
      real ( kind = 8 ) xyz_0(3),r0

      real *8, allocatable :: errs(:),cms(:,:),rads(:)
      real *8 thet,phi
      complex * 16 zpars(3)
      integer numit,niter
      character *200 title,fname,fname1,fname2,fname3,fname_base,fname_out

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima,zk,zvec(3)
      complex *16 alpha_rhs
      character (len=10), dimension(2) :: scalar_names

      integer count1

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)


      done = 1
      pi = atan(done)*4

      zk=1.0d0

      zpars(1) = zk 
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

      xyz_in(1) = 0.11d0
      xyz_in(2) = 0.0d-5
      xyz_in(3) = 0.37d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 20.1d0

      norder = 16
      iref = 2
      fname_base = '../geometries_go3/cow_new'
      write(fname,'(a,a,i2.2,a,i1,a)') trim(fname_base),'_o', &
        norder,'_r0',iref,'.go3'
      write(fname_out,'(a,a,i2.2,a,i1,a)') trim(fname_base),'_o', &
        norder,'_r0',iref,'_quad.go3'
      write(fname1,'(a,i2.2,a,i1,a)') '../vtk_files/cow_new_o', &
        norder,'_r0',iref,'.vtk'
      print *, "fname=",trim(fname)
            
      call open_gov3_geometry_mem(fname,npatches,npts)

      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,&
         iptype,npts,srcvals,srccoefs,wts)

      allocate(cms(3,npatches),rads(npatches))


      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,cms,rads)

      call prin2('srcvals=*',srcvals(3,npts-12:npts),12)
      call prin2('normals=*',srcvals(12,npts-12:npts),12)

      call prin2('srcvals=*',srcvals(3,1:12),12)
      call prin2('normals=*',srcvals(12,1:12),12)

      xmin = minval(srcvals(1,:))
      xmax = maxval(srcvals(1,:))
      
      ymin = minval(srcvals(2,:))
      ymax = maxval(srcvals(2,:))
      
      zmin = minval(srcvals(3,:))
      zmax = maxval(srcvals(3,:))

      
      
      dlam = 2*pi/real(zk)

      dx = (xmax-xmin)/dlam
      dy = (ymax-ymin)/dlam
      dz = (zmax-zmin)/dlam

      rmin = minval(rads)
      rmax = maxval(rads)

      print *, "dx=",dx
      print *, "dy=",dy
      print *, "dz=",dz
      
      print *, "rmin=",rmin
      print *, "rmax=",rmax

      print *, "rrat=",rmax/rmin
      

      call plot_surface_info_all(dlam,npatches,norders,ixyzs,iptype, &
        npts,srccoefs,srcvals,trim(fname1),'a')
      
      call trimesh_to_quadmesh(fname_base,norder,iref,fname_out)


      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars,& 
     &srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri,& 
     &triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,& 
     &ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,&
     &npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 0
        vmax = 2*pi
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),&
     &nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,& 
     &iptr3,iptr4, norder,'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,&
     &npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end


      subroutine test_exterior_pt(npatches,norder,npts,srcvals,&
     &srccoefs,wts,xyzout,isout)
!
!
!  this subroutine tests whether the pt xyzin, is
!  in the exterior of a surface, and also estimates the error
!  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
!  centered at the interior point. Whether a point 
!  is in the interior or not is tested using Gauss' 
!  identity for the flux due to a point charge
!
!
!  input:
!    npatches - integer
!       number of patches
!    norder - integer
!       order of discretization
!    npts - integer
!       total number of discretization points on the surface
!    srccoefs - real *8 (9,npts)
!       koornwinder expansion coefficients of geometry info
!    xyzout -  real *8 (3)
!       point to be tested
!
!  output: 
!    isout - boolean
!      whether the target is in the interior or not
!

      implicit none
      integer npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
      real *8 ra,ds
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2


      zk = 0

      ra = 0



      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          call h3d_sprime(xyzout,srcvals(1,i),dpars,zk,ipars,val)
          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-3) isout = .false.
      if(abs(ra).le.1.0d-3) isout = .true.

      return
      end

   

subroutine record_solution_nrccie(a_vect,ns,norder,filename)
implicit none

    !List of calling arguments
    integer, intent(in) :: ns, norder
    complex ( kind = 8 ), intent(in) :: a_vect(2*ns)
    character (len=*) filename

    !List of local variables
    integer ( kind = 4 ) umio,count1,count2,flag,npols
    integer :: ierror,ntri
    complex ( kind = 8 ) ima
		
    ima=(0.0d0,1.0d0)
    npols = (norder+1)*(norder+2)/2
    ntri=ns/npols

    open(8, FILE=filename,STATUS='REPLACE')
    write(8,*) norder
    write(8,*) ntri
    do count1=1,3*ns
        write(8,*) real(a_vect(count1))
    enddo
    do count1=1,3*ns
        write(8,*) imag(a_vect(count1))
    enddo

    close (8)

return
end subroutine record_solution_nrccie



subroutine 	get_rhs_em_nrccie_pec_sphere(P0,vf,alpha,ns,srcvals,zk,RHS)
implicit none
!
!  This function obtains the right hand side for the NRCCIE formulation
!  for the integral boundary equation:
!
!    J/2-M_{k}[J]+alpha·nxnx(ikS_{k}[J]-gradS_{k}[rho]) =
!      = nxH_inc + alpha nxnxE_inc
!    rho/2+S'_{k}[rho]-ikS_{k}[J]+alpha(divS_{k}[J]-ikS_{k}[rho]) = 
!      = n·E_inc
!
!  input:
!    P0 - real * 8 (3)
!      location of the source point at the exterior region
!      WARNING! notice that this formulation uses a representation theorem
!      for the incoming field in the interior region (MFIE) therefore
!      therefore it only works for incoming fields generated by sources in
!      the exterior region (or at infinity like plane waves)
!
!    vf - complex *16(3)
!      Orientation of the magnetic and electric dipoles located at P0 
!
!    alpha - complex *16
!      parameter in the combined formulation
!   
!    ns - integer
!      total number of points on the surface
!
!    srcvals - real *8(12,ns)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    zk - complex *16
!      Helmholtz parameter 
!
!  output:
!    RHS - complex  *16(3*ns)
!      right hand side
!      RHS(1:ns) - first component of  nxH_inc + alpha nxnxE_inc along
!       the srcvals(4:6,i) direction
!      RHS(ns+1:2*ns) - second component of  nxH_inc + alpha nxnxE_inc
!       along the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      RHS(2*ns+1:3*ns) - normal component of the electric field n·E_inc
!

	!List of calling arguments
	integer ( kind = 4 ), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk,alpha
	complex ( kind = 8 ), intent(out) :: RHS(3*ns)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
	integer count1
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
  complex ( kind = 8 ) ima
		
  data ima/(0.0d0,1.0d0)/

	allocate(E(3,ns), H(3,ns))

!	call fieldsED(zk,P0,srcvals,ns,E,H,vf,0)
!	call fieldsMD(zk,P0,srcvals,ns,E,H,vf,1)

  do count1=1,ns
    E(1,count1)=exp(ima*zk*srcvals(3,count1))
    E(2,count1)=0.0d0
    E(3,count1)=0.0d0

    H(1,count1)=0.0d0
    H(2,count1)=exp(ima*zk*srcvals(3,count1))
    H(3,count1)=0.0d0

!    write (*,*) srcvals(1,count1)*srcvals(4,count1)+srcvals(2,count1)*srcvals(5,count1)+srcvals(3,count1)*srcvals(6,count1)
  enddo

	do count1=1,ns	
      call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)
      RHS(count1)=-DOT_PRODUCT(rv,H(:,count1))+&
	     &alpha*DOT_PRODUCT(ru,E(:,count1))
      RHS(ns+count1)=DOT_PRODUCT(ru,H(:,count1))+&
       &alpha*DOT_PRODUCT(rv,E(:,count1))
      RHS(2*ns+count1)=DOT_PRODUCT(srcvals(10:12,count1),E(:,count1))
	enddo

return
end subroutine get_rhs_em_nrccie_pec_sphere







subroutine 	get_rhs_em_nrccie_pec_sphere2(P0,vf,alpha,ns,srcvals,zk,RHS)
implicit none
!
!  This function obtains the right hand side for the NRCCIE formulation
!  for the integral boundary equation:
!
!    J/2-M_{k}[J]+alpha·nxnx(ikS_{k}[J]-gradS_{k}[rho]) =
!      = nxH_inc + alpha nxnxE_inc
!    rho/2+S'_{k}[rho]-ikS_{k}[J]+alpha(divS_{k}[J]-ikS_{k}[rho]) = 
!      = n·E_inc
!
!  input:
!    P0 - real * 8 (3)
!      location of the source point at the exterior region
!      WARNING! notice that this formulation uses a representation theorem
!      for the incoming field in the interior region (MFIE) therefore
!      therefore it only works for incoming fields generated by sources in
!      the exterior region (or at infinity like plane waves)
!
!    vf - complex *16(3)
!      Orientation of the magnetic and electric dipoles located at P0 
!
!    alpha - complex *16
!      parameter in the combined formulation
!   
!    ns - integer
!      total number of points on the surface
!
!    srcvals - real *8(12,ns)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    zk - complex *16
!      Helmholtz parameter 
!
!  output:
!    RHS - complex  *16(3*ns)
!      right hand side
!      RHS(1:ns) - first component of  nxH_inc + alpha nxnxE_inc along
!       the srcvals(4:6,i) direction
!      RHS(ns+1:2*ns) - second component of  nxH_inc + alpha nxnxE_inc
!       along the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      RHS(2*ns+1:3*ns) - normal component of the electric field n·E_inc
!

	!List of calling arguments
	integer ( kind = 4 ), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk,alpha
	complex ( kind = 8 ), intent(out) :: RHS(3*ns)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
	integer count1
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
  complex ( kind = 8 ) ima
		
  data ima/(0.0d0,1.0d0)/

	allocate(E(3,ns), H(3,ns))

!	call fieldsED(zk,P0,srcvals,ns,E,H,vf,0)
!	call fieldsMD(zk,P0,srcvals,ns,E,H,vf,1)

  do count1=1,ns
    E(1,count1)=exp(ima*zk*srcvals(2,count1))
    E(2,count1)=0.0d0
    E(3,count1)=0.0d0

    H(1,count1)=0.0d0
    H(2,count1)=0.0d0
    H(3,count1)=-exp(ima*zk*srcvals(2,count1))

!    write (*,*) srcvals(1,count1)*srcvals(4,count1)+srcvals(2,count1)*srcvals(5,count1)+srcvals(3,count1)*srcvals(6,count1)
  enddo

	do count1=1,ns	
      call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)
      RHS(count1)=-DOT_PRODUCT(rv,H(:,count1))+&
	     &alpha*DOT_PRODUCT(ru,E(:,count1))
      RHS(ns+count1)=DOT_PRODUCT(ru,H(:,count1))+&
       &alpha*DOT_PRODUCT(rv,E(:,count1))
      RHS(2*ns+count1)=DOT_PRODUCT(srcvals(10:12,count1),E(:,count1))
	enddo

return
end subroutine get_rhs_em_nrccie_pec_sphere2





subroutine 	get_rhs_em_nrccie_pec_sphere3(P0,vf,alpha,ns,srcvals,zk,RHS)
implicit none
!
!  This function obtains the right hand side for the NRCCIE formulation
!  for the integral boundary equation:
!
!    J/2-M_{k}[J]+alpha·nxnx(ikS_{k}[J]-gradS_{k}[rho]) =
!      = nxH_inc + alpha nxnxE_inc
!    rho/2+S'_{k}[rho]-ikS_{k}[J]+alpha(divS_{k}[J]-ikS_{k}[rho]) = 
!      = n·E_inc
!
!  input:
!    P0 - real * 8 (3)
!      location of the source point at the exterior region
!      WARNING! notice that this formulation uses a representation theorem
!      for the incoming field in the interior region (MFIE) therefore
!      therefore it only works for incoming fields generated by sources in
!      the exterior region (or at infinity like plane waves)
!
!    vf - complex *16(3)
!      Orientation of the magnetic and electric dipoles located at P0 
!
!    alpha - complex *16
!      parameter in the combined formulation
!   
!    ns - integer
!      total number of points on the surface
!
!    srcvals - real *8(12,ns)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    zk - complex *16
!      Helmholtz parameter 
!
!  output:
!    RHS - complex  *16(3*ns)
!      right hand side
!      RHS(1:ns) - first component of  nxH_inc + alpha nxnxE_inc along
!       the srcvals(4:6,i) direction
!      RHS(ns+1:2*ns) - second component of  nxH_inc + alpha nxnxE_inc
!       along the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      RHS(2*ns+1:3*ns) - normal component of the electric field n·E_inc
!

	!List of calling arguments
	integer ( kind = 4 ), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk,alpha
	complex ( kind = 8 ), intent(out) :: RHS(3*ns)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
	integer count1
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
  complex ( kind = 8 ) ima
		
  data ima/(0.0d0,1.0d0)/

	allocate(E(3,ns), H(3,ns))

!	call fieldsED(zk,P0,srcvals,ns,E,H,vf,0)
!	call fieldsMD(zk,P0,srcvals,ns,E,H,vf,1)

  do count1=1,ns
    E(1,count1)=0.0d0
    E(2,count1)=0.0d0
    E(3,count1)=exp(-ima*zk*srcvals(1,count1))

    H(1,count1)=0.0d0
    H(2,count1)=exp(-ima*zk*srcvals(1,count1))
    H(3,count1)=0.0d0

!    write (*,*) srcvals(1,count1)*srcvals(4,count1)+srcvals(2,count1)*srcvals(5,count1)+srcvals(3,count1)*srcvals(6,count1)
  enddo

	do count1=1,ns	
      call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)
      RHS(count1)=-DOT_PRODUCT(rv,H(:,count1))+&
	     &alpha*DOT_PRODUCT(ru,E(:,count1))
      RHS(ns+count1)=DOT_PRODUCT(ru,H(:,count1))+&
       &alpha*DOT_PRODUCT(rv,E(:,count1))
      RHS(2*ns+count1)=DOT_PRODUCT(srcvals(10:12,count1),E(:,count1))
	enddo

return
end subroutine get_rhs_em_nrccie_pec_sphere3


