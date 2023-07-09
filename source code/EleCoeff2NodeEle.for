C ======================================================================
c
c User Subroutine EleCoeff2NodeEle for sbfem
c All rights of reproduction or distribution in any form are reserved.
c By Yang Yang (PhD), yanghhu@foxmail.com
c Parameter:
c          Coord: the coordinates of node | coord(dim,nnode)
c          kcond    : the permeability
c          kss    : the specific storage coefficient
c Return:
c          the element coefficient matrices - e0, e1, e2, m0
c
C ======================================================================
      subroutine EleCoeff2NodeEle(coord,kcondx,kcondy,kss,E0,E1,E2,M0)
      implicit double precision(A-H,O-Z)
      double precision coord(2,2),C1(2,1),C2(2,1),K(2,2)
      double precision Q0(1,1),Q1(1,1),Q2(1,1),C1T(1,2),C2T(1,2)
      double precision Q0_temp(1,2),Q1_temp(1,2),Q2_temp(1,2)
      double precision E0(2,2),E1(2,2),E2(2,2),M0(2,2)
      double precision meanX,meanY,deltaX,deltaY,area,Jb
      double precision kcondx,kss,kcondy
c
c     Initiation Q0,Q1,Q2
c
      Q0(1,1)=0.0D0
c
      Q1(1,1)=0.0D0
c
      Q2(1,1)=0.0D0
c  
c     Initiation Q0_temp,Q1_temp,Q2_temp
      do i=1,2
        Q0_temp(1,i)=0.0D0
      end do
c
      do i=1,2
        Q1_temp(1,i)=0.0D0
      end do
c
      do i=1,2
        Q2_temp(1,i)=0.0D0
      end do
c
c     Initiation e0,e1,e2
      do i=1,2
         do j=1,2
            e0(i,j)=0.0D0
         end do
      end do
c   
      do i=1,2
         do j=1,2
            e1(i,j)=0.0D0
         end do
      end do
c
      do i=1,2
         do j=1,2
            e2(i,j)=0.0D0
         end do
      end do
c
      deltaX=coord(1,2)-coord(1,1)
      deltaY=coord(2,2)-coord(2,1)
      Jb=0.5d0*(coord(1,1)*coord(2,2)-coord(1,2)*coord(2,1))
      meanX=0.5d0*(coord(1,1)+coord(1,2))
      meanY=0.5d0*(coord(2,1)+coord(2,2))
      area=coord(1,1)*coord(2,2)-coord(1,2)*coord(2,1)
c
      if (area.lt.1e-10) then 
         write(*,*) "negative area"
      end if
c
      C1(1,1)=deltaY*0.5d0
      C1(2,1)=-deltaX*0.5d0
      C1T=transpose(C1)
c
      C2(1,1)=-meanY
      C2(2,1)=meanX
      C2T=transpose(C2)
c
c     Seepage matrix K
      K(1,1)=kcondx
      K(1,2)=0.0d0
      K(2,1)=0.0d0
      K(2,2)=kcondy

c     Calculate Q0,Q1,Q2
      Q0_temp=0.5d0*(1.0d0/Jb)*matmul(C1T,K)
      Q0=matmul(Q0_temp,C1)
c
      Q1_temp=0.5d0*(1.0d0/Jb)*matmul(C2T,K)
      Q1=matmul(Q1_temp,C1)
c     
      Q2_temp=0.5d0*(1.0d0/Jb)*matmul(C2T,K)
      Q2=matmul(Q2_temp,C2)
c
c     Calculate E0,E1,E2
      e0(1:1,1:1)=2.0d0*Q0
      e0(1:1,2:2)=Q0
      e0(2:2,1:1)=Q0
      e0(2:2,2:2)=2.0d0*Q0
      e0=(2.0d0/3.0d0)*e0
c
      e1(1:1,1:1)=-1.0d0*Q1-(1.0d0/3.0d0)*Q0
      e1(1:1,2:2)=-1.0d0*Q1+(1.0d0/3.0d0)*Q0
      e1(2:2,1:1)=Q1+(1.0d0/3.0d0)*Q0
      e1(2:2,2:2)=Q1-(1.0d0/3.0d0)*Q0
c
      e2(1:1,1:1)=Q2+(1.0d0/3.0d0)*Q0
      e2(1:1,2:2)=-1.0d0*Q2-(1.0d0/3.0d0)*Q0
      e2(2:2,1:1)=-1.0d0*Q2-(1.0d0/3.0d0)*Q0
      e2(2:2,2:2)=Q2+(1.0d0/3.0d0)*Q0
c
      M0(1,1)=2.0d0
      M0(1,2)=1.0d0
      M0(2,1)=1.0d0
      M0(2,2)=2.0d0
c
      M0=M0*(kss*Jb/3.0d0)
      return
      end
