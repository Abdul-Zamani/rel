!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This program computes relativistic shifts to
!atomic core ionization potential using the 
!BSE. 
!
!Some references
!https://aip.scitation.org/doi/pdf/10.1063/1.1733236 Silverman/Scherr
!https://journals.aps.org/pr/pdf/10.1103/PhysRev.112.1649 Perekis
!https://journals.aps.org/pr/pdf/10.1103/PhysRev.139.A619 Bagus
!https://iopscience.iop.org/article/10.1088/0022-3700/18/5/008/pdf Schirmer&Barth
!https://link.springer.com/book/10.1007/978-3-662-12869-5 (Bethe/Salpeter)
!
! AZ 2/23/2022
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program rel
      implicit none
      integer :: i,j
      real :: z
      real :: a0r1,a1r1,a2r1,a3r1
      real :: a0r12,a1r12,a2r12,a3r12
      real :: a0p4,a1p4,a2p4,a3p4
      real :: r1,r12,p4,e,e2,fine,au2ev
      real :: a0e2,a1e2,a2e2,a3e2
      real, dimension(20) :: eVec,e2Vec

      a0r1 = 1.0
      a1r1 = -0.667
      a2r1 = 0.174
      a3r1 = 0.010

      a0r12 = 0.125
      a1r12 = -0.2430
      a2r12 = 0.186
      a3r12 = -0.073

      a0p4 = 5
      a1p4 = -4.0917
      a2p4 = 1.96
      a3p4 = -0.63

      a0e2 = 0.16282
      a1e2 = -0.2400
      a2e2 = 0.1179
      a3e2 = -0.022



      fine =  0.0072973525693
      au2ev = 27.211386245988

      r1=0
      r12=0
      p4=0
      e=0

      j=0
      do i=1,20
        z=i
        r1=(z**3.0)*(a0r1*z**(0.0) + a1r1*z**(-1.0) &
          + a2r1*z**(-2.0) + a3r1*z**(-3.0))
        r12=(z**3.0)*(a0r12*z**(0.0) + a1r12*z**(-1.0) &
          + a2r12*z**(-2.0) + a3r12*z**(-3.0))
        p4=(z**4.0)*(a0p4*z**(0.0) + a1p4*z**(-1.0) &
          + a2p4*z**(-2.0) + a3p4*z**(-3.0))

        e2 = (z**2.0)*(a0e2*z**(0.0) + a1e2*z**(-1.0) &
          + a2e2*z**(-2.0) + a3e2*z**(-3.0))

        e= (fine**2.0)*((((-1.0)*z**4.0)/8.0) + (p4/4.0) - &
           (z*r1) - (r12))



        j=j+1
        !print*,'z',z,'j',j
        !print*,'r1',r1,'r12',r12,'p4',p4,'e2',e2
        eVec(j) = e
        e2Vec(j) = (fine**2.0)*((e2/2.0))

      enddo 

      write(*,*)'Relativistic Shifts to IP (Z=1,20) w/ and w/o E2'
      do i=1,20
        print*,'Z=',i,'Ej=',eVec(i)*au2ev,'Ej(E2)=',&
        (eVec(i)+e2Vec(i))*au2ev
      enddo


      end program rel
