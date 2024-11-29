      program ljmd
c***********************************************************************
c     
c     program je napravljem prema ideji 
C
C     CCP5 Summer School Lennard Jones Molecular Dynamics Program
c     
c     copyright Daresbury Laboratory
c     author W. Smith 2003
c     
c***********************************************************************

      implicit none

      integer matms,natms,ngrid
      parameter (matms=3,ngrid=128,natms=4*matms*matms*matms)

      logical ltraj,lpos
      integer i,j,k,l,n,istep,nstep,neql,nav,isamp,nrep
      real*8 x(natms),y(natms),z(natms)
      real*8 vx(natms),vy(natms),vz(natms)
      real*8 fx(natms),fy(natms),fz(natms)
      real*8 gor(ngrid)
      real*8 sid,mx,my,mz,scl,tst,rct,bol,time,pi,rrr,nrm
      real*8 tot,totav,pot,vir,kin,tmp,tmpav,den,tem,dtm,vol
      real*8 prs,prsav,box,dr
      real*8 xn(4),yn(4),zn(4)
      data xn/-.25d0,0.25d0,-.25d0,0.25d0/
      data yn/-.25d0,0.25d0,0.25d0,-.25d0/
      data zn/-.25d0,-.25d0,0.25d0,0.25d0/
      data pi/3.1415926536d0/


c     ulazni parametri simulacije 

      read(*,*)nstep
      read(*,*)neql
      read(*,*)isamp
      read(*,*)den
      read(*,*)tem
      read(*,*)dtm
      read(*,*)ltraj
      read(*,*)lpos
      write(*,*)'Broj vremenskih koraka        ',nstep
      write(*,*)'Broj ekvilibracijskih koraka  ',neql
      write(*,*)'Interval sakupljanja          ',isamp
      write(*,*)'Broj atoma                    ',natms
      write(*,*)'Gustoća sustava               ',den
      write(*,*)'Temperatura sustava           ',tem
      write(*,*)'Duljina vremenskog koraka     ',dtm
      write(*,*)'Parametar (T) - ispiši trajektorije u HISTORY.',ltraj
c      write(*,*)'Parametar (T) upisuje položaje, (F) brzine u HISTORY. ',lpos

c     početak simulacije - inicijalizacija

      nav=0
      tot=0.d0
      totav=0.d0
      prsav=0.d0
      tmpav=0.d0
      bol=1.5d0
      tst=dtm*0.5d0
      vol=dble(natms)/den
      box=vol**(1.d0/3.d0)
      sid=box/dble(matms)
      rct=0.5d0*box
      dr=rct/dble(ngrid)
      nrep=max(1,nstep/20)

c     početni položaji atoma - fcc rešetka 

      k=1
      do n=1,4
        do i=1,matms
          do j=1,matms
            do l=1,matms
              x(k)=(dble(l)-0.5d0+xn(n))*sid-0.5d0*box
              y(k)=(dble(j)-0.5d0+yn(n))*sid-0.5d0*box
              z(k)=(dble(i)-0.5d0+zn(n))*sid-0.5d0*box
              k=k+1
            enddo
          enddo
        enddo
      enddo

c     inicijalizacija polja za računanje rdf

      do i=1,ngrid
        gor(i)=0.d0
      enddo

c     početne brzine

      mx=0.d0
      my=0.d0
      mz=0.d0
      call random_number(vx)
      call random_number(vy)
      call random_number(vz)
      do k=1,natms
        vx(k)=vx(k)-0.5d0
        vy(k)=vy(k)-0.5d0
        vz(k)=vz(k)-0.5d0
        mx=mx+vx(k)
        my=my+vy(k)
        mz=mz+vz(k)
      enddo
      mx=mx/dble(natms)
      my=my/dble(natms)
      mz=mz/dble(natms)

c     skaliranje da ukupna količina gibanja bude 0

      kin=0.d0
      do k=1,natms
        vx(k)=vx(k)-mx
        vy(k)=vy(k)-my
        vz(k)=vz(k)-mz
        kin=kin+vx(k)*vx(k)+vy(k)*vy(k)+vz(k)*vz(k)
      enddo

c     skaliranje za zadanu temperaturu 

      kin=0.5d0*kin
      tmp=kin/(bol*dble(natms))
      scl=sqrt(tem/tmp)
      do k=1,natms
        vx(k)=scl*vx(k)
        vy(k)=scl*vy(k)
        vz(k)=scl*vz(k)
      enddo

c     otvaramo datoteku za ispis statistike termodinamičkih veličina

      open(9,file="STA")

c     otvaramo datoteku u koju se upisuju konfiguracije

      if(ltraj)open(11,file="TRJ")

c     računanje početnih sila --- funkcija forces - ova funkcija računa uz sile i energije, tlak i RDF

      call forces
     x  (natms,pot,vir,rct,box,dr,x,y,z,fx,fy,fz,gor)

c    integracijski algoritam -  velocity verlet algorithm
C    računanje jednadžba gibanja 

      do istep=1,nstep

        pot=0.d0
        kin=0.d0
        vir=0.d0
        
c     prvi dio 
        
        do k=1,natms
          vx(k)=vx(k)+tst*fx(k)
          vy(k)=vy(k)+tst*fy(k)
          vz(k)=vz(k)+tst*fz(k)
          x(k)=x(k)+dtm*vx(k)
          y(k)=y(k)+dtm*vy(k)
          z(k)=z(k)+dtm*vz(k)
        enddo
        
c     računanje sila 
        
        call forces
     x    (natms,pot,vir,rct,box,dr,x,y,z,fx,fy,fz,gor)
        
c    međukorak za računanje jednadžbi gibanja
        
        do k=1,natms
          vx(k)=vx(k)+tst*fx(k)
          vy(k)=vy(k)+tst*fy(k)
          vz(k)=vz(k)+tst*fz(k)
          kin=kin+vx(k)*vx(k)+vy(k)*vy(k)+vz(k)*vz(k)
        enddo
        kin=0.5d0*kin

c     periodični granični uvjeti 
        
        do k=1,natms
          x(k)=x(k)-box*nint(x(k)/box)
          y(k)=y(k)-box*nint(y(k)/box)
          z(k)=z(k)-box*nint(z(k)/box)
        enddo
        
c     računanje termodinamički veličina sustava 

        tmp=kin/(bol*natms)
        tot=kin+pot
        prs=(2.d0*kin-vir)/(3.d0*vol)
        if(istep.gt.neql)then
          nav=nav+1
          tmpav=((tmpav*(nav-1))/nav)+(tmp/nav)
          totav=((totav*(nav-1))/nav)+(tot/nav)
          prsav=((prsav*(nav-1))/nav)+(prs/nav)
        endif
        
c     checkpoint ispis --- dobro je uvijek staviti neki ispis --- kontrolira izvršavanje programa !!!!!!!

        time=dtm*istep
        if(mod(istep,nrep).eq.0)
     x    write(*,'(1x,1p,4e16.8)')time,tot,tmp,prs
        if(istep.eq.neql)write(*,*)'Kraj ekvilibracije!!!'

c     upis u datoteku za statistiku - STA datoteka

        write(9,'(1p,4e16.8)')time,tot,tmp,prs

        if (mod(istep,isamp).eq.0)then
          
c     upis konfiguracija - TRJ datoteka

          if(ltraj)then
            write(11,'(i10,f12.4,f12.7)')natms,time,box
            if(lpos)then
              do i=1,natms
                write(11,'(3f12.7)')x(i),y(i),z(i)
              enddo
            else
              do i=1,natms
                write(11,'(3f12.7)')vx(i),vy(i),vz(i)
              enddo
            endif
          endif

c     skaliranje za zadanu temperaturu --- najosnovniji termostat !!!!!!
        
          if(istep.le.neql)then
            scl=sqrt(tem/tmp)
            do k=1,natms
              vx(k)=scl*vx(k)
              vy(k)=scl*vy(k)
              vz(k)=scl*vz(k)
            enddo
          endif

        endif

c     odbaciti dio koji ulazi u ekvilibraciju za računanje rdf 

        if(istep.le.neql)then
          do i=1,ngrid
            gor(i)=0.d0
          enddo
        endif

      enddo

c   konačne srednje veličine energije, temperature i tlaka 

      write(*,*)'ukupna srednja energija, temperatura i tlak:'
      write(*,'(1x,i10,1p,3e16.8)')nav,totav,tmpav,prsav

c ispisuje zadnju konfiguraciju 
      open(8,file="XYZ")
      write(8,'(i10,f12.7)')natms,box
      write(8,'(a)')'Lennard Jonesium'
      do i=1,natms
        write(8,'(a4,3f12.7)')'Ar  ',x(i),y(i),z(i)
      enddo

c     računanje rdf 

      open(10,file="RDF")
      nrm=1.d0/(2.d0*pi*dr*den*dble(nstep-neql)*dble(natms-1))
      do i=1,ngrid
        rrr=(dble(i)-0.5d0)*dr
        gor(i)=nrm*gor(i)/(rrr**2+dr**2/12.d0)
        write(10,'(2f12.7)')rrr,gor(i)
      enddo
      if(ltraj)close(11)      
      close(8)
      close(9)
      close(10)
      write(*,*)'Job done'

      end

      subroutine forces
     x  (natms,pot,vir,rct,box,dr,x,y,z,fx,fy,fz,gor)
c***********************************************************************
c     
c     program za računanje sile za lennard-jones potencijal    
c     
c***********************************************************************

      implicit none

      logical lrdf
      integer natms,i,j,ix
      real*8 pot,vir,rct,box,rc2,alp,bet,dx,dy,dz,dr
      real*8 sg6,sg12,gam,rrs,rr0,rr3,rrr,rsq

      real*8 x(*),y(*),z(*),fx(*),fy(*),fz(*),gor(*)

      lrdf=.true.
      pot=0.d0
      vir=0.d0
      rc2=rct*rct
      alp=24.d0*(2.d0*rct**(-12)-rct**(-6))/rct
      bet=4.d0*(rct**(-12)-rct**(-6))+alp*rct
      
      do i=1,natms
        fx(i)=0.d0
        fy(i)=0.d0
        fz(i)=0.d0
      enddo

      do i=1,natms-1
        do j=i+1,natms

          dx=x(j)-x(i)
          dy=y(j)-y(i)
          dz=z(j)-z(i)
c     mimimalna udaljenost - periodični uvjeti 
          dx=dx-box*nint(dx/box)
          dy=dy-box*nint(dy/box)
          dz=dz-box*nint(dz/box)
c     udaljenost čestica
          rsq=dx*dx+dy*dy+dz*dz
c     sfereni r_cutoff - definira radijus do kojeg uključujemo čestice u dvočestičnu interakciju 
          if (rsq.lt.rc2)then
            rrr=sqrt(rsq)
c     rdf  račun 
            if(lrdf)then
              ix = int(rrr/dr)+1
              gor(ix)=gor(ix)+1.d0
            endif
c     računanje potencijalne energije 
            rrs=1.d0/rsq
            rr0=rrr*rrs
            rr3=rr0*rrs
            sg6=rr3*rr3
            sg12=sg6*sg6
            pot=pot+4.d0*(sg12-sg6)+alp*rrr-bet
c     računanje sile 
            gam=24.d0*rrs*(2.d0*sg12-sg6)-alp*rr0
            if(gam.gt.1.d6)gam=1.d6
            vir=vir-gam*rsq
c projekcije na sve osi (x,y,z) i uključeno da na i-tu česticu djeluje sila od j-te, a na j-tu ista sila samo suprotnog predznaka 
            fx(i)=fx(i)-gam*dx
            fy(i)=fy(i)-gam*dy
            fz(i)=fz(i)-gam*dz
            fx(j)=fx(j)+gam*dx
            fy(j)=fy(j)+gam*dy
            fz(j)=fz(j)+gam*dz

          endif
        enddo
      enddo
      return
      end
