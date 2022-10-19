c=======================================================================
c Ray tracing code using explicit formula.                             c
c By Yan-Rong Li                                                       c
c liyanrong@ihp.ac.cn                                                  c
c 2008-04-21 to 2008-04-30                                             c                                            
c=======================================================================

c***********************************************************************
c to cal the radii in the disk where the photon is emitted.
c input:
c      astar:
c      theta:
c      alpha:
c      beta:
c
c output:
c      rem:
c      lambda:
c      Q:
c      sr:    
c***********************************************************************

      subroutine raytrac(astar,theta,alpha,beta,rem,lambda,Q,sr)
      implicit none
      include 'const.f'
      real*8 ellf
      
      real*8 astar,theta,theta_obs      
      real*8 alpha,beta,rem,lambda,Q,srr
      integer sr
      
      real*8 aspin,Rr,rootmax
      real*8 mup2,mum2,mmu,cn,dn,sn,taumu,cmu,amu,cnobs,sinobs,p
      
      real*8 eacc
      integer i,signb
      eacc=1.0d-4
      
      theta_obs=theta*pi/180.0d0
      aspin=astar/2.0d0

      if(beta.gt.0.0d0)then
      signb=1
      elseif(beta.lt.0.0d0)then
      signb=-1
      else
      signb=0
      endif
      
      lambda=-alpha*dsin(theta_obs)
      Q=beta**2.0d0+(alpha**2.0d0-aspin**2.0d0)*(dcos(theta_obs))**2.0d0
c      write(*,*)"======================================================"
c      write(*,*)"Lambda=",lambda,",  Q=",Q      
      if((Q.lt.0.0d0).or.((Q.eq.0.0d0)
     &.and.(abs(lambda).le.abs(aspin))))then
c      write(*,*)"No solution!"
      rem=0.0d0
      sr=-1.0d0
      return
      endif
      
      if(aspin.ne.0.0d0)then

      mup2=0.5d0/aspin**2.0d0
     &    *(dsqrt((lambda**2.0d0+Q-aspin**2.0d0)**2.0
     &      +4.0d0*aspin**2.0d0*Q)
     &    -(lambda**2.0d0+Q-aspin**2.0d0))
     
c      mum2=0.5d0/aspin**2.0d0
c     &    *(dsqrt((lambda**2.0d0+Q-aspin**2.0d0)**2.0
c     &      +4.0d0*aspin**2.0d0*Q)
c     &    +(lambda**2.0d0+Q-aspin**2.0d0))
     
      mum2=Q/aspin/aspin/mup2
           
      mmu=mup2/(mup2+mum2)
      mmu=dsqrt(mmu) 
      amu=dsqrt(aspin**2.0d0*(mup2+mum2))
    
c      write(*,*)"======================================================"      
c      write(*,*)"mu+=",dsqrt(mup2),", mu-=",dsqrt(mum2)
c     &          ,", mu_obs=",dcos(theta_obs)
      cnobs=dcos(theta_obs)/dsqrt(mup2)
      cnobs=min(cnobs,1.0d0)
      cnobs=dacos(cnobs)
      taumu=(ellf(pi/2.0,mmu)+signb*ellf(cnobs,mmu))/amu            
      else

      cmu=dsqrt(Q/(lambda**2.0+Q))
c      write(*,*)"======================================================"
c      write(*,*)"mu+=",cmu,", mu_obs=",dcos(theta_obs)
      sinobs=min(dcos(theta_obs)/cmu,1.0)
      taumu=pi/2.0d0+signb*(pi/2.0d0-dasin(sinobs))
      taumu=taumu/dsqrt(lambda**2.0+Q)  
      end if  
c note that conversion from r_g to r_s. 
c      call taumucal(p,astar,alpha*2.0,beta*2.0,theta_obs)
c      write(*,*)taumu,p*2.0
c      taumu=p*2.0d0 
    
      call remsolve(taumu/2.0d0,rem,srr,astar,lambda*2.0d0,Q*4.0d0)   
      sr=int(srr) 
      rem=rem/2.0d0
      end subroutine raytrac
c*********************************************************
c subroutine taumucal()
c provided by Prof. Yuan
c*********************************************************      
      subroutine taumucal(taumu,astar,alpha,beta,theta_obs)
      implicit none
      real*8 taumu,astar
      real*8 pmax1,pmax2,a,lambda,q,theta_obs,dth_em,
     *     xmiu_em,xmiu1,xmiu2,xmiu1_sq,xmiu2_sq,xmth,psi_obs,
     *     xmiu12,xmiu_obs,thetain,ql,PI,q_min,p0,EPS,err1,err2,
     *     alpha,beta,theta_em,dth_obs,psi_em,erro,rh
      real*8 ellf,ellk
      parameter(EPS=1.0d-5)      
      
      PI=dacos(0.0d0)*2.0
      rh=1.0+dsqrt(1.-astar*astar)
      theta_em=PI/2.0d0
      
      a=astar
      
      taumu=0.0d0
      
      lambda=-alpha*dsin(theta_obs)
      Q=beta**2.-(a*dcos(theta_obs))**2.
     *     +(alpha*dcos(theta_obs))**2.

      xmiu_em=0.0d0
      xmiu_obs=dcos(theta_obs)
c...
c... Case I: a=0
c...
      if(dabs(a).eq.0.0d0) then
      
      if(lambda.eq.0.0d0.and.q.eq.0.0d0) then
       taumu=1./rh
      endif
      p0=dsqrt(q+lambda**2.0d0)
      xmiu1=dsqrt(q/(q+lambda*lambda))
      if(xmiu1.gt.1.0d0) xmiu1=1.0d0
      if(xmiu_em.gt.xmiu1) then
      write(*,*)"xmiu_em.gt.xmiu1"
      
      return
      endif
      
      if(xmiu_obs.gt.xmiu1) then      
      write(*,*)"xmiu_obs.gt.xmiu1"     
      return
      endif
      
      dth_obs=xmiu_obs/xmiu1
      if(dth_obs.ge.1.0d0) dth_obs=1.0d0
      psi_obs=dacos(dth_obs)
      
      dth_em=xmiu_em/xmiu1
      if(dth_em.ge.1.0d0) dth_em=1.0d0
      psi_em=dacos(dth_em)
      
      if(beta.gt.0.0d0) then
      taumu=(psi_obs+psi_em)/p0
      else if(beta.lt.0.0d0)then
      taumu=(psi_em-psi_obs)/p0
      else
      taumu=psi_em/p0
      endif
      
      return
      endif

c...
c... Case II: a\=0
c...
c...

      if(dabs(lambda).eq.0.0d0) then
      xmiu1_sq=1.0d0
      xmiu1=1.0d0
      else
      xmiu1_sq=(dsqrt((lambda**2.+q-a**2.)**2.+4.0*a*a*q)-
     *        (lambda**2.+q-a**2.))/(2.0*a*a)
      xmiu1=dsqrt(xmiu1_sq)
      if(xmiu1.ge.1.0d0) xmiu1=1.0d0
      endif

      erro=dabs((xmiu_em-xmiu1)/xmiu_em)
      if(erro.lt.1.0d-3.and.xmiu_em.gt.xmiu1) xmiu1=xmiu_em

      if(q.ge.0.0d0) then
      if(beta.eq.0.0d0) xmiu1=xmiu_obs
      if(xmiu_em.gt.xmiu1) return
      if(xmiu_obs.gt.xmiu1) return

      xmiu2_sq=q/a**2./xmiu1_sq
      xmiu2=dsqrt(xmiu2_sq)
      xmth=xmiu1/dsqrt(xmiu1_sq+xmiu2_sq)
      p0=dsqrt(a*a*(xmiu1_sq+xmiu2_sq))
      dth_em=xmiu_em/xmiu1
      dth_obs=xmiu_obs/xmiu1
      psi_em=ellf(dacos(dth_em),xmth)
      psi_obs=ellf(dacos(dth_obs),xmth)
      
      if(beta.gt.0.0d0) then
      taumu=(psi_obs+psi_em)/p0
      else if(beta.lt.0.0d0)then
      taumu=(psi_em-psi_obs)/p0
      else
      taumu=psi_em/p0
      endif
      return
      
      endif

      if(q.lt.0.0d0) then
c      write(*,*)"No solution!"
      return
      endif      
      end subroutine taumucal
c********************************************************
c subroutine remsolve(p,rem,sr,a,lambda,q)
c provided by Prof. Yuan Ye-Fei
c note that the units is r_g=0.5r_s
c********************************************************      
      subroutine remsolve(p,rem,sr,a,lambda,q)
      implicit none
      real*8 a,lambda,q,theta_jet,p,rem,sr,aa,
     *    phi_rh,ksi_rh,p_rh,phi_infty,ksi_infty,
     *    p_infty,xa,xb,p0,s_rh,s_rem
      real*8 b,c,d,e,f,d1,d2,d3,d4,xtheta,ra,rb,
     *    rc,rd,r1,r2,r3,r4
      real*8 xm4,uu,emmc,sn,cn,dn,u,v,w
      real*8 xm2,omegar,eta0,c1,c2,dd1,dd2
      real*8 cbrt,ellf,ellf2
      real*8 z1,z2,rms,dp,rh
      real*8 xc,xd,ld,l1,l2,p1,p2,pp,q1,q2,qq,fr,r

      if(p.le.0.0d0) then
      rem=1.0d30
      sr=1.0d0
      return
      endif

      if(a.eq.0.0d0.and.lambda.eq.0.0d0.and.q.eq.0.0d0) then
      rem=1/p
      sr=1.0d0
      return
      endif

      rh=1.0+dsqrt(1.-a**2)

      c=(a-lambda)**2+q
      d=2.0/3.0*(q+lambda**2-a**2)
      e=9.0/4.0*d**2-12.0*a**2*q
      f=-27.0/4.0*d**3-108.0*a**2*q*d+108.0*c**2
	
      d1=f**2-4.0*e**3
      if(d1.ge.0.0) then
      aa=1./3.*(cbrt(0.5*(f-dsqrt(d1)))+cbrt(0.5*(f+dsqrt(d1))))
      else
      d2=dabs(d1)
      d3=dsqrt(f**2+d2)
      xtheta=dacos(f/d3)
      aa=1./3.*cbrt(0.5*d3)*2.0*dcos(xtheta/3.0)
      endif		

      b=dsqrt(aa+d)

c-----------------------------------------------------------------
      dd1=-aa+2.0*d-4.0*c/b
      dd2=-aa+2.0*d+4.0*c/b

      u=0.5*b
      w=0.5*dsqrt(dabs(dd1))
      v=0.5*dsqrt(dabs(dd2))

c... 
c... Case I: four real roots (dd1>0 & dd2>0)
c...
      if(dd1.gt.0.0d0.and.dd2.gt.0.0d0) then 
      ra=u-w
      rb=u+w
      rc=-u-v
      rd=-u+v
      call sort_root(ra,rb,rc,rd,r1,r2,r3,r4)

      if(r2.eq.r1) pause'r2=r1'

      xm4=dsqrt(((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3)))
      p0=2.0/dsqrt((r2-r4)*(r1-r3))
	
      phi_infty=dasin(dsqrt((r2-r4)/(r1-r4)))
      ksi_infty=ellf(phi_infty,xm4)
      p_infty=p0*ksi_infty

      emmc=1.0-xm4*xm4

c... r1 >= rh
      if(r1.ge.rh) then  ! r1 > rh
      dp=2.0*p_infty
c... rem         
      if(p.gt.dp) then
      rem=0.0d0
      sr=1.0d0
      else
      uu=(p_infty-p)/p0
      call sncndn(uu,emmc,sn,cn,dn)
      rem=(r1*(r2-r4)-r2*(r1-r4)*sn**2)/(r2-r4-(r1-r4)*sn**2)
      if(p.lt.p_infty) then
      sr=1.0d0
      else
      sr=-1.0d0
      endif
      endif    
c... r1 < rh
      else
      phi_rh=dasin(dsqrt((r2-r4)*(rh-r1)/(r1-r4)/(rh-r2)))
      ksi_rh=ellf(phi_rh,xm4)
      p_rh=p0*ksi_rh

      dp=p_infty-p_rh
      if(p.gt.dp) then
      rem=0.0d0
      sr=1.0d0
      else
      uu=(p_infty-p)/p0
      call sncndn(uu,emmc,sn,cn,dn)
      rem=(r1*(r2-r4)-r2*(r1-r4)*sn**2)/
     *    (r2-r4-(r1-r4)*sn**2)
      sr=1.0d0
      endif
      endif
      endif

c...
c... Case II: two real roots and two complex roots (Ref. Li etal 2005)
c...

      if(dd1.lt.0.0d0.and.dd2.gt.0.0d0) then  ! two real roots
      r4=-u-v                  
      r3=-u+v         
      xa=dsqrt((r3-u)**2.+w*w)   ! w = v in Li etal 2005
      xb=dsqrt((r4-u)**2.+w*w) 
      xm2=dsqrt(((xa+xb)**2.-(r3-r4)**2.)/(4.*xa*xb))
      
      p0=1./dsqrt(xa*xb)
      
      phi_infty=dacos((xa-xb)/(xa+xb))
      
      ksi_infty=ellf2(phi_infty,xm2)
      p_infty=p0*ksi_infty
      
      emmc=1.0-xm2*xm2

c... r3 >= rh

      if(r3.ge.rh) then  ! r3 > rh
      
      dp=2.0*p_infty
c... rem 

      if(p.gt.dp) then
      rem=0.0d0
      sr=1.0d0
      
      else
      
      uu=(p_infty-p)/p0
      call sncndn(uu,emmc,sn,cn,dn)
      rem=(r4*xa-r3*xb-(r4*xa+r3*xb)*cn)/
     &    ((xa-xb)-(xa+xb)*cn)
     
      if(p.lt.p_infty) then
      sr=1.0d0
      else
      sr=-1.0d0
      endif
      
      endif
c... r3 < rh
      else
      
      phi_rh=dacos(((xa-xb)*rh+r3*xb-r4*xa)/
     &        ((xa+xb)*rh-r3*xb-r4*xa))
      ksi_rh=ellf2(phi_rh,xm2)
      p_rh=p0*ksi_rh
      
      dp=p_infty-p_rh
      if(p.gt.dp) then
      rem=0.0d0
      sr=1.0d0
      else
      
      uu=(p_infty-p)/p0
      call sncndn(uu,emmc,sn,cn,dn)
      rem=(r4*xa-r3*xb-(r4*xa+r3*xb)*cn)/
     &     ((xa-xb)-(xa+xb)*cn)
      sr=1.0d0
      endif
      endif
      endif
c...

c...  
c... Case III: four complex roots
c...

      if(dd1.lt.0.0d0.and.dd2.lt.0.0) then
      xb=2.0*u
      xc=u*u+w*w
      xd=u*u+v*v
      ld=dsqrt((xd-xc)**2+2.0*xb**2*(xd+xc))
      l1=0.5*((xc-xd)-ld)/xb
      l2=0.5*((xc-xd)+ld)/xb
      p1=l1**2-xb*l1+xc
      p2=l1**2+xb*l1+xd
      q1=l2**2-xb*l2+xc
      q2=l2**2+xb*l2+xd
      pp=max(q1/p1,q2/p2)
      qq=min(q1/p1,q2/p2)
      xm2=dsqrt((pp-qq)/pp)
      emmc=1.0-xm2*xm2
      
      p0=(l2-l1)/dsqrt(pp*p1*p2)
      phi_infty=dasin(dsqrt(1./(1.+qq)))
      ksi_infty=ellf(phi_infty,xm2)
      p_infty=p0*ksi_infty

      s_rh=(rh-l2)/(rh-l1)
      phi_rh=dasin(s_rh*dsqrt(1./(s_rh**2+qq)))
      ksi_rh=ellf(phi_rh,xm2)
      p_rh=p0*ksi_rh
      if(p.le.(p_infty-p_rh)) then
      uu=(p_infty-p)/p0
      call sncndn(uu,emmc,sn,cn,dn)
      s_rem=dsqrt(qq)*sn/dsqrt(1.-sn**2)
      rem=(s_rem*l1-l2)/(s_rem-1.)
      sr=1.0d0
      else
      rem=0.0d0
      sr=1.0d0
      endif
      endif
      
      if(dd1.eq.0.0d0.or.dd2.eq.0.0d0) then
      print*,'equal real roots'
      pause
      endif
      
      return 
      end
      
      subroutine sort_root(ra,rb,rc,rd,r1,r2,r3,r4)
      implicit none
      real*8 ra,rb,rc,rd,r1,r2,r3,r4,rr(4),r
      integer n,i,j
      
      rr(1)=ra
      rr(2)=rb
      rr(3)=rc
      rr(4)=rd
      
      n=4
      
      do j=2,n
      r=rr(j)
      do i=j-1,1,-1
      if(rr(i).le.r) goto 10
      rr(i+1)=rr(i)
      enddo
      i=0
10    rr(i+1)=r

      enddo

      r1=rr(4)
      r2=rr(3)
      r3=rr(2)
      r4=rr(1)
      return 
      end
	
      double precision function cbrt(x)
      implicit none
      double precision x
      if(x.eq.0.0d0) then 
      cbrt=0.0d0
      else 
      cbrt=x/dabs(x)*dabs(x)**(1./3.)
      endif
      return
      end      
c*****************************************
c special functions from Numerical Recipes.
c*****************************************     
      SUBROUTINE laguer(a,m,x,its)
      INTEGER m,its,MAXIT,MR,MT
      REAL*8 EPSS
      COMPLEX*8 a(m+1),x
      PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)
      INTEGER iter,j
      REAL*8 abx,abp,abm,err,frac(MR)
      COMPLEX*8 dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=cmplx(0.,0.)
        f=cmplx(0.,0.)
        abx=abs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=abs(b)+abx*err
11      continue
        err=EPSS*err
        if(abs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.*f/b
          sq=sqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          abp=abs(gp)
          abm=abs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.) then
            dx=m/gp
          else
            dx=exp(cmplx(log(1.+abx),float(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
12    continue
      pause 'too many iterations in laguer'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.   
      FUNCTION ellf(phi,ak)
      REAL*8 ellf,ak,phi
CU    USES rf
      REAL*8 s,rf
      s=dsin(phi)
      ellf=s*rf(dcos(phi)**2,(1.-s*ak)*(1.+s*ak),1.0d0)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.

      FUNCTION ellf2(phi,ak)
      implicit none
      double precision ellf2,ak,phi,PI,ellk
CU    USES rf
      double precision s,rf
      PI=dacos(0.0d0)*2.0
      s=dsin(phi)

      ellf2=s*rf(dcos(phi)**2,(1.-s*ak)*(1.+s*ak),1.d0)
      if(phi.gt.PI/2.) ellf2=2*ellk(ak)-ellf2
      return
      END

      FUNCTION ellk(ak)
      implicit none
      double precision ellk,ak
CU    USES rf
      double precision s,rf
      s=1.0d0
      ellk=s*rf(0.0d0,(1.-s*ak)*(1.+s*ak),1.d0)
      return
      END      

      FUNCTION rf(x,y,z)
      REAL*8 rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08,TINY=1.5e-38,BIG=3.E37,THIRD=1./3.,
     *C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
      REAL*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=dsqrt(xt)
        sqrty=dsqrt(yt)
        sqrtz=dsqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/dsqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#   
	
      SUBROUTINE sncndn(uu,emmc,sn,cn,dn)
      REAL*8 cn,dn,emmc,sn,uu,CA
      PARAMETER (CA=0.0003d0)
      INTEGER i,ii,l
      REAL*8 a,b,c,d,emc,u,em(13),en(13)
      LOGICAL bo
      emc=emmc
      u=uu
      if(emc.ne.0.)then
        bo=(emc.lt.0.)
        if(bo)then
          d=1.-emc
          emc=-emc/d
          d=dsqrt(d)
          u=d*u
        endif
        a=1.
        dn=1.
        do 11 i=1,13
          l=i
          em(i)=a
          emc=dsqrt(emc)
          en(i)=emc
          c=0.5*(a+emc)
          if(abs(a-emc).le.CA*a)goto 1
          emc=a*emc
          a=c
11      continue
1       u=c*u
        sn=dsin(u)
        cn=dcos(u)
        if(sn.eq.0.)goto 2
        a=cn/sn
        c=a*c
        do 12 ii=l,1,-1
          b=em(ii)
          a=c*a
          c=dn*c
          dn=(en(ii)+a)/(b+a)
          a=c/b
12      continue
        a=1./dsqrt(c**2+1.)
        if(sn.lt.0.)then
          sn=-a
        else
          sn=a
        endif
        cn=c*sn
2       if(bo)then
          a=dn
          dn=cn
          cn=a
          sn=sn/d
        endif
      else
        cn=1./dcosh(u)
        dn=cn
        sn=dtanh(u)
      endif
      return
      END

