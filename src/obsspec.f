c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                 supported by Prof. Wang Jian-Min
c                          and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-01-
c======================================================================= 

c To cal the observed spectrum with GR effects included.

c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<c
c                                                                      c
c compile this code, following files are needed:                       c
c                                                                      c
c raytrac.f   ######   for ray tracing                                 c
c obsvar.f    ######   for variables definition                        c
c const.f     ######   for some constants                              c
c                                                                      c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>c

c
c Thu, Aug 28,2008
c rewrite the interpolation. Thanks Prof.Yuan's help.

      include 'raytrac.f'
      
      program obsspec
      implicit none
      include 'const.f'
      include 'obsvar.f'
      real*8 y_at_x
      integer ngli,nglo,nphi
      parameter(ngli=100,nglo=100,nphi=200)  
      real*8 alpha,beta,rem,theta_obs,lambda,q
     &,     mesh_rho(ngli+nglo),mesh_phi(nphi)
     &,     mesh_rem(ngli+nglo,nphi,4)
     &,     mesh_drho(ngli+nglo)
      namelist/obs/theta_obs
      real*8 rho,phi,rli1,rli2,rlo1,rlo2,drli,drlo,dphi
      real*8 ds,gr,gp,omg,vel,log_rem

      real*8 redshift,nulnu(nnu),nu_unred(nnu),nu_obs(nnu)
     &,      flux_obs(nnu),dnu,nu_min,nu_max
      real*8 aspin,ametric,delta,vr,rs,rmin,rmax,m,zmax,zmin
      integer sr
      
      integer i,j,kk,ird,np

      real*8 mdot,gammaa,visf
      namelist/adaf/m,mdot,alpha,beta,gammaa,astar,visf
      open(unit=8,file='./data/datain.txt',status='old')
      read(8,adaf)
      write(*,adaf)
      read(8,obs)
      write(*,obs)
      close(8)
      pause
      open(unit=10,file='./data/specobs.dat',status='unknown')

      aspin=astar/2.0d0
      rmin=(1.0d0+sqrt(1.0d0-astar*astar))/2.0d0+0.2d0 
      rmax=1.0d3 ! 1000Rs          
      rs=2.95d5*m
      
      call read_sol()
      call read_flux()
                   
      rli1=rmin
      rli2=10.0d0

      rlo1=rli2
      rlo2=rmax

c from rli1 to rli2, we use the linear mesh grid

c      drli=(rli2-rli1)/(ngli)
      drli=dlog10(rli2/rli1)/(ngli)
      drlo=dlog10(rlo2/rlo1)/nglo

      dphi=(2.0*pi)/nphi
      
c divide the mesh (alpha,beta)

      do i=1,ngli
c      mesh_rho(i)=rli1+(i-0.5d0)*drli
c      mesh_drho(i)=drli
      mesh_rho(i)=rli1*10.0**((i-0.5d0)*drli)
      mesh_drho(i)=rli1*10.0**(i*drli)*(1.0d0-10.0**(-drli))
      enddo
      
      do i=1,nglo
      mesh_rho(i+ngli)=rlo1*10.0**((i-0.5d0)*drlo)
      mesh_drho(i+ngli)=rlo1*10.0**(i*drlo)*(1.0d0-10.0**(-drlo))
      enddo
      do i=1,nphi
      mesh_phi(i)=(i-0.5d0)*dphi
      enddo
    
      do i=1,nglo+ngli
        do j=1,nphi
         alpha=mesh_rho(i)*dcos(mesh_phi(j))
         beta=mesh_rho(i)*dsin(mesh_phi(j))
         call raytrac(astar,theta_obs,alpha,beta,rem,lambda,Q,sr)
         mesh_rem(i,j,1)=rem
         mesh_rem(i,j,2)=lambda
         mesh_rem(i,j,3)=q
         mesh_rem(i,j,4)=sr
c         write(20,'(5f15.7,I5)')mesh_rho(i),rem,alpha,beta
        enddo
      enddo

      nu_max=22.0d0
      nu_min=9.0d0
      dnu=(nu_max-nu_min)/nnu
      do kk=1,nnu
      nu_obs(kk)=nu_min+(kk-1)*dnu
c      nu_obs(kk)=log_nu(kk)
      enddo

c begin integrate over the coordinates (alpha,beta)
      do kk=1,nnu
      nulnu(kk)=0.0d0+1.0d-150
      enddo
      
      zmax=-1.0d10
      zmin=1.0d10
      
      do 40 i=1,nglo+ngli
      
c      write(*,*)"i=",i

      do 30 j=1,nphi
      
c      if(i.eq.4.and.j.eq.85)then
c      write(*,*)i,j,nu_obs(2),nulnu(2),flux_obs(2)
c      endif
      
      rem=mesh_rem(i,j,1)
      if(rem.lt.rmin)goto 30
      if(rem.gt.rmax)goto 30
      
      lambda=mesh_rem(i,j,2)
      q=mesh_rem(i,j,3)
      sr=mesh_rem(i,j,4)
      
      log_rem=dlog10(rem)
      
c interprate      
      vel=y_at_x(log_rem,log_rd,log_v,nrd,np)
      vel=10.0d0**(vel)
      gp=y_at_x(log_rem,log_rd,log_gamp,nrd,np)
      gp=10.0d0**(gp)

c we did not use log10(omega) since omega will be negative when a<0.      
      omg=y_at_x(log_rem,log_rd,omega,nrd,np)
      
      if(vel.ge.1.0d-6)then
      gr=1.0d0/sqrt(1.0d0-vel*vel)
      else
      gr=1.0d0+0.5d0*vel*vel
      endif
      if(gr.lt.1.0d0) gr=1.0d0
      if(gp.lt.1.0d0) gp=1.0d0      
      
       
      ds=mesh_drho(i)*mesh_rho(i)*dphi  ! the area element.
      
      ametric=1.0d0+aspin**2.0d0/rem**2.0d0 + aspin**2.0d0/rem**3.0d0
      delta=1.0d0-1.0d0/rem + aspin**2.0d0/rem**2.0d0
       
      vr=((rem**2.0d0+aspin**2.0d0)-aspin*lambda)**2.0d0
     &    -rem**2.0d0*delta*((lambda-aspin)**2.0d0+q)

c       vr=rem**4.0d0+(aspin*aspin-lambda**2.0-Q)*rem*rem
c     &   +(Q+(lambda-aspin)**2.0d0)*rem-aspin**2.0d0*Q
     
      if(vr.lt.0.0d0)then
      write(*,*)rem,mesh_rho(i)*dcos(mesh_phi(j)),
     &          mesh_rho(i)*dsin(mesh_phi(j))
      write(*,*)lambda,q,sr,vr
      write(*,*)aspin,delta
      vr=0.0d0
      pause
      endif
      
      Vel=-vel
      vr=sr*dsqrt(vr)
      redshift=gr*gp*ametric**0.5d0/delta**0.5d0
     & *(1.0d0-Vel*vr/gp/rem**2.0d0/ametric**0.5d0
     & -omg*lambda)-1.0d0  
      if(redshift.lt.-1.0d0)pause
      
      if(zmin.gt.redshift)zmin=redshift
      if(zmax.lt.redshift)zmax=redshift
      
c      redshift=0.0d0
      do kk=1,nnu
      nu_unred(kk)=nu_obs(kk)+dlog10(1.0d0+redshift)
      enddo
      
      call interp_nu(log_rem,nu_unred,flux_obs)

      do kk=1,nnu
      nulnu(kk)=nulnu(kk)+10.0d0**(nu_obs(kk))*flux_obs(kk)
     &                 *ds/(1.0d0+redshift)**3.0d0*4.0*PI
      enddo
      
c      write(*,*)i,j,log10(nulnu(180))
c      if(i.eq.15.and.j.eq.50)then
c      write(*,*)i,j,rem,nu_obs(nnu),nu_unred(nnu),flux_obs(nnu)
c      write(*,*)i,j,nu_obs(2),nu_unred(2),flux_obs(2)
c      endif

30    continue
40    continue

      
      do kk=1,nnu
      write(10,'(e15.7,a,e15.7)')nu_obs(kk),"  "
     &,          dlog10(nulnu(kk))+2.0d0*dlog10(rs)
      write(*,*)nu_obs(kk),dlog10(nulnu(kk))+2.0d0*dlog10(rs)   
      enddo
      write(*,*)zmin,zmax
      close(10)
      end program obsspec
      
c********************************************************
c subroutine interp_nu(log_rem,nu_unred,flux_obs)
c********************************************************     
      subroutine interp_nu(log_rem,nu_unred,flux_obs)
      implicit none
      include 'obsvar.f'
      real*8 y_at_x
      real*8 nu_unred(nnu),flux_obs(nnu),log_rem
      real*8 rd_int(nrd),flux_nrd(nrd),log_flux_int(nnu)
      real*8 nus,fnu
      integer k,m,inu,nrdp,nnup,n_nearest_pt
      real*8 error

c interpolte the r direction firstly.      
      do k=1,nnu
      
      do m=1,nrd
      flux_nrd(m)=log_flux(m,k)
      enddo
      log_flux_int(k)=y_at_x(log_rem,log_rd,flux_nrd,nrd,n_nearest_pt)
c      write(*,*)n_nearest_pt,10.0**(log_rem)
      enddo

c then the frequency direction
c if the frequency exceed the bounday of log_nu, there will appear errors
c in calling the function y_at_x. We just set its value to be the boundary.
      
      do 110 k=1,nnu
      
      nus=nu_unred(k)
      
      if(nus.lt.log_nu(1))then
      fnu=log_flux_int(1)
      goto 111
      endif
      if(nus.gt.log_nu(nnu))then
      fnu=log_flux_int(nnu)
      goto 111
      endif          
      fnu=y_at_x(nus,log_nu,log_flux_int,nnu,n_nearest_pt)
      
111   if(fnu.lt.-150.0d0)fnu=-150.0d0
      flux_obs(k)=10.0d0**(fnu)
    
110   continue   
      return
      end subroutine interp_nu
c********************************************************
c subroutine read_flux()
c********************************************************
      subroutine read_flux()
      implicit none
      include 'const.f'
      include 'obsvar.f'
      integer ird,i,j,stat,inu

      do ird=1,nrd
      
      open(unit=40,file="./data/spec/spec"//char(ird/100+48)//
     &char(mod(ird,100)/10+48)//char(mod(ird,10)+48)//".txt",
     &form="formatted",status='old')
     
      do inu=1,nnu
      read(40,*,iostat=stat)log_nu(inu),log_flux(ird,inu)
      if(stat.gt.0)goto 90 
      nu(inu)=10.0**(log_nu(inu))
c
c note the factor 2.0*PI
      log_flux(ird,inu)=log_flux(ird,inu)
      flux(ird,inu)=10.0**(log_flux(ird,inu))
      enddo
      
90    close(40)

      do j=inu,nnu
      log_flux(ird,j)=-150.0d0
      flux(ird,j)=1.0d-150
      enddo

      enddo
      return
      end subroutine read_flux
c********************************************************
c subroutine read_sol()
c********************************************************
      subroutine read_sol()
      implicit none
      include 'obsvar.f'
      integer i
      real*8 rmin
      rmin=(1.0d0+sqrt(1.0d0-astar*astar))/2.0d0+0.2d0
      open(unit=20,file='./data/sol_for_spec.dat',status='old')
      do i=1,nobs
      read(20,*,end=50)rd(i),gamr(i),gamp(i),omega(i),v(i)
      if(rd(i).lt.rmin)goto 50
      log_rd(i)=dlog10(rd(i))
      log_gamr(i)=dlog10(gamr(i))
      log_gamp(i)=dlog10(gamp(i))
      omega(i)=omega(i)
      log_v(i)=dlog10(v(i))
      
      enddo
50    close(20)
      nrd=i-1
      write(*,*)"NRD=",nrd
      return
      end subroutine read_sol
c********************************************************
c  function y_at_x()
c********************************************************
      function y_at_x(logx,log_x_tab,log_y_tab,n_tab,n_nearest_pt)
      implicit none
      integer n_tab,n_nearest_pt
      double precision y_at_x,log_x_tab(n_tab),log_y_tab(n_tab),
     *                   logx,logy,interp
      logy=interp(logx,log_x_tab,log_y_tab,n_tab,n_nearest_pt)
      y_at_x=logy
      return 
      end 
c********************************************************
c SUBROUTINE hunt(xx,n,x,jlo)
c********************************************************  
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.

c/*************************************************************************/
c/* Driver for the interpolation routine. First we find the tab. point    */
c/* nearest to xb, then we interpolate using four points around xb.       */
c/*************************************************************************/
      function interp(xb,xp,yp,np,n_nearest_pt)
      implicit none
      integer  np,n_nearest_pt,k,m
c..                    k--index of 1st point
c..                    m--degree of interpolation
      double precision xp(np),yp(np),xb,interp,DBL_EPSILON
      parameter(DBL_EPSILON=1.0d-15)

      m=4

      call hunt(xp,np,xb,n_nearest_pt)

c      if(n_nearest_pt.eq.0)then
c      interp=yp(1)
c      return
c      endif
c      if(n_nearest_pt.eq.np)then
c      interp=yp(np)
c      return
c      endif  
          
      k=min(max((n_nearest_pt)-(m-1)/2,1),np+1-m)

      if(xb.eq.xp(k).or.xb.eq.xp(k+1).or.xb.eq.xp(k+2).or.xb.eq.xp(k+3))
     *      xb = xb + DBL_EPSILON

      interp = (xb-xp(k+1))*(xb-xp(k+2))*(xb-xp(k+3))*yp(k)/
     1        ((xp(k)-xp(k+1))*(xp(k)-xp(k+2))*(xp(k)-xp(k+3)))
     2        +(xb-xp(k))*(xb-xp(k+2))*(xb-xp(k+3))*yp(k+1)/
     3        ((xp(k+1)-xp(k))*(xp(k+1)-xp(k+2))*(xp(k+1)-xp(k+3)))
     4        +(xb-xp(k))*(xb-xp(k+1))*(xb-xp(k+3))*yp(k+2)/
     5        ((xp(k+2)-xp(k))*(xp(k+2)-xp(k+1))*(xp(k+2)-xp(k+3)))
     6        +(xb-xp(k))*(xb-xp(k+1))*(xb-xp(k+2))*yp(k+3)/
     7        ((xp(k+3)-xp(k))*(xp(k+3)-xp(k+1))*(xp(k+3)-xp(k+2)))
      return
      end
	
