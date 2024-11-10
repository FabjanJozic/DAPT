*********************************************************
* FCC CRYSTAL CONFIG  (SEPT 92)
*********************************************************
      program start
      implicit real*8 (a-h,o-z)
      parameter (ndim=10976, xfcc=0.1547)
      dimension rx(ndim),ry(ndim),rz(ndim)
      dimension vx(ndim),vy(ndim),vz(ndim)

      write(*,*)'INPUT RHO* , N  (NPARTICLE=4*N^3  :N=3->NP=108 etc..)'
      read(*,*)rho,n 
       np=4*n**3
       write(*,*)'NUMBER OF PARTICLE=',NP   
        
      

*find crystal parameters 
       vol=(np/rho)
       xl=vol**0.33333333
       diam=1./xl

       
       write(*,'(''Total volume             '',f12.5)')vol
       write(*,'(''Total number of particle '',i5)')np
       write(*,'(''Distance per side        '',f12.5)')xl
       write(*,'(''Particle diameter(red)   '',f12.5)')diam

C
C     SET UP STARTING FCC LATTICE (NSP MUST BE OF THE FORM 4*N**3)
      IN=(NP/4)**(1./3.)+0.1
      IN2=2*IN
      UNIT=1.0/FLOAT(IN)
      SHIFT=0.2
C
c      write(*,*) IN,IN2,UNIT
      IL=1
      XA=-1.0-SHIFT
      DO 20 JX=1,IN2
      XA=XA+UNIT
      YA=-1.0-SHIFT
      DO 20 JY=1,IN2
      YA=YA+UNIT
      ZA=-1.0-SHIFT
      DO 20 JZ=1,IN2
      ZA=ZA+UNIT
      IXYZ=JX+JY+JZ
      IF((IXYZ/2)*2.NE.IXYZ) GO TO 20
      RX(IL)=XA*0.5
      RY(IL)=YA*0.5
      RZ(IL)=ZA*0.5
      IL=IL+1
   20 CONTINUE

        
       IPART=IL-1
 10    continue
       miss=np-ipart
       if(miss.ne.0) write(*,'(//,
     &'' ATTENTION WITH ACTUAL NPART'',2i4,//)')ipart,np

*test overlap
       ibug=0
       diam2=diam*diam
       do 100 i=1,ipart-1
       do 100 j=i+1,ipart
       if(i.ne.j) then
         rxij=rx(i)-rx(j)
         ryij=ry(i)-ry(j)
         rzij=rz(i)-rz(j)
              RXIJ = RXIJ - ANINT ( RXIJ )
              RYIJ = RYIJ - ANINT ( RYIJ )
              RZIJ = RZIJ - ANINT ( RZIJ )
         r=RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ  
         if(r.lt.diam2) then
          ibug=ibug+1
          zoz=sqrt(r)/diam
          write(*,'(''OVERLAP BETWEEN '',2i4,2f10.5,f15.9)')
     &i,j,sqrt(r),zoz,abs(zoz-1)
          write(*,'(3f10.4)')rx(i),ry(i),rz(i)
          write(*,'(3f10.4)')rx(j),ry(j),rz(j)
c          write(*,'(3f10.4)')rx(i)/diam,ry(i)/diam,rz(i)/diam
c          write(*,'(3f10.4)')rx(j)/diam,ry(j)/diam,rz(j)/diam
          endif
         
       endif
 100   continue

       write(*,*)'TOTAL NUMBER OF OVERLAPPING PARTICLES : ',ibug

*assign random opposite velocities to pair of particles
       k=1
       vmax=0.1
  111  continue
       k=k+1
       xk=real(k)
       i1=ipart*ranf(xk)
       k=k+1
       yk=real(k)
       i2=ipart*ranf(yk)
       if(i1.ne.0.and.i2.ne.0.and.i1.ne.i2) then
        vx(i1)=vmax*ranf(rx(i1))
        vy(i1)=vmax*ranf(ry(i1))
        vz(i1)=vmax*ranf(rz(i1))
        vx(i2)=-vx(i1)   !! opposite to insure zero momentum
        vy(i2)=-vy(i1)   !! opposite to insure zero momentum
        vz(i2)=-vz(i1)   !! opposite to insure zero momentum
       endif
       if(k.lt.ipart) goto 111

*displacement parameters
       drmax=0.15*diam

       open(1,file='mc.xyz')
       write(1,*)ipart , drmax, rho
       do 50 i=1,ipart
 50    write(1,'(6f15.6)')rx(i),ry(i),rz(i),vx(i),vy(i),vz(i)
** 50    write(1,*)rx(i),ry(i),rz(i)

       close(1)
       end
*******************************************************************
** RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.         **
**                                                               **
**                 ***************                               **
**                 **  WARNING  **                               **
**                 ***************                               **
**                                                               **
** GOOD RANDOM NUMBER GENERATORS ARE MACHINE SPECIFIC.           **
** PLEASE USE THE ONE RECOMMENDED FOR YOUR MACHINE.              **
*******************************************************************
        FUNCTION RANF ( DUMMY )
        implicit real*8(a-h,o-z)

        INTEGER     L, C, M
        PARAMETER ( L = 1029, C = 221591, M = 1048576 )

        INTEGER     SEED
c       REAL        DUMMY
        SAVE        SEED
        DATA        SEED / 0 /


        SEED = MOD ( SEED * L + C, M )
        RANF = REAL ( SEED ) / M

        RETURN
        END


