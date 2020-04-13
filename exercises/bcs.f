**********************************************************************
      Subroutine BCS(Np,Gp,E,delta,lambda,Epc,vk2)
**********************************************************************
****** Newtons method for solving a pair of non-linear equations *****
      Implicit Double Precision (A-H,O-Z)
      Dimension E(364),vk2(364)
      Real*8 Jac, lambda
      x=0.86D0
      y0=E(Np/2)
      y=y0

5     Do 10 i=1,100
      iFlag=0
      Call Funcs(Np,Gp,x,y,F,dFx,dFy,G,dGx,dGy,E,iFlag)
      If(iFlag.EQ.0) Goto 15
        Dx=F*dGy-dFy*G
        Dy=dFx*G-F*dGx
        Jac=dFx*dGy-dFy*dGx
        x1=x-Dx/Jac
        y1=y-Dy/Jac
        If((DAbs(Dx/Jac).LT.0.00001).AND.(DAbs(Dy/Jac).LT.0.001))Then
          delta=DAbs(x1)
          lambda=y1
        Goto 20
        Endif
        x=x1
        y=y1
10    Continue

15    Write(*,*)'BCS equations are not converging, Delta set to 0.0'
      delta=0.0D0
      lambda=y0
      Epc=0.0D0
      Return
20    Continue
      Sekvk2=0.0D0
      Sek=0.0D0
      Svk4=0.0D0
      S1=0.0D0
      Do 30 i = 1,Np
      Dr =((E(i)-lambda)**2+delta**2)**0.5D0
        vk2(i)=0.5D0*(1.0D0 - ((E(i)-lambda)/Dr))
        Sekvk2=Sekvk2+ E(i)*vk2(i)
        Svk4=Svk4+(vk2(i))**2
        If (i.GT.(Np/2)) Goto 30
          Sek=Sek+ E(i)
          S1=S1+ 1.0D0
30    Continue
      Epc = 2.0D0*(Sekvk2-Sek) - delta**2/Gp - Gp*(Svk4-S1)
      End

***********************************************************************
      Subroutine Funcs(Np,Gp,x,y,F,dFx,dFy,G,dGx,dGy,E,iFlag)
      Implicit Double Precision (A-H,O-Z)
      Dimension E(364)
      iFlag=0
      F=(-1.0D0)*DFloat(Np)
      G=(-1.0D0)*2.0D0/Gp
      dFx=0.0D0
      dFy=0.0D0
      dGx=0.0D0
      dGy=0.0D0
      Do 10 I=1,Np
        Dr =((E(I)-y)**2+x**2)**0.5D0
        If (DAbs(Dr).LT.0.000001) Return
        If (DAbs(Dr).GT.100000) Return
        dDr=Dr*((E(I)-y)**2+x**2)
        F=F+(1.0D0 - (E(I)-y)/Dr)
        dFx=dFx + x*(E(I)-y)/dDr
        dFy=dFy + x**2/dDr
        G=G+1.0D0/Dr
        dGx=dGx - x/dDr
        dGy=dFy + (E(I)-y)/dDr
10    Continue
      iFlag=1
      End
