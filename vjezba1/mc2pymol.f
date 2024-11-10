c This codes transforms my MC files into a Pymol compatible .xyz
c format file, with atoms drawn and carbon atoms; and cell size
c scaled appropriately
      program vmdrender
      character*70 nom
     
      cdiam=4  ! approx size of carbon atom in Angstrom
      shift=0  !!0.5 ! shift all coord by half box size
      write(*,*)'Input the MC file name to render in Pymol'
      read(*,'(a)')nom
      open(1,file=nom)
       read(1,*)n,dr,rho
       n1=n
*compute scaling factor
      cfact=cdiam*(n1/rho)**(1./3.)
      write(*,'("Lbox=",f12.5)')cfact
c      m=n1  !-1   ! for some reason ????
c      if(ione.ne.1) m=m+n2   !! 2018
      open (2,file="dataMC.xyz")
      write(2,'(i10,"  LBOX=",f10.5)')n1,cfact
      do i=1,n1
       read(1,*)x,y,z
       x=(x+shift)*cfact
       y=(y+shift)*cfact
       z=(z+shift)*cfact
       write(2,'("C     ",3f12.7)')x,y,z
      enddo
      close(1)
      close(2)
      write(*,*)'VMD file written in "dataMC.xyz" '
      end
