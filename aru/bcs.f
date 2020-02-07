      Subroutine BCS(Np,Gp,x0,y0,E,Epc,delta)
*****************Newtons method for******************
*******solving a pair of non-linear equations********
*****************************************************
      Dimension E(364),vk2(364)
      Real Jac, lambda
      x=x0
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
        If((Abs(Dx/Jac).LT.0.00001).AND.(Abs(Dy/Jac).LT.0.001))Then
          delta=Abs(x1)
          lambda=y1
          Goto 20
        Endif
        x=x1
        y=y1
10    Continue
15    Write(*,*)'BCS equations are not converging, Epc set to 0.0'
      Epc=0.0
      Return
20    Continue
      Sekvk2=0.0
      Sek=0.0
      Svk4=0.0
      S1=0.0
      Do 30 i = 1,Np
        Dr =((E(i)-lambda)**2+delta**2)**0.5
        vk2(i)=0.5*(1.0 - ((E(i)-lambda)/Dr))
        Sekvk2=Sekvk2+ E(i)*vk2(i)
        Svk4=Svk4+(vk2(i))**2
        If (i.GT.(Np/2)) Goto 30
          Sek=Sek+ E(i)
          S1=S1+ 1.0
30    Continue
      write(*,*)Sekvk2,Svk4,Sek,S1
      Epc = 2.0*(Sekvk2-Sek) - delta**2/Gp - Gp*(Svk4-S1)
      Write(*,'(3(5X,A,F10.5))')
     &      'Delta =',delta,'Lambda=',lambda,'Epc =',Epc

      End

      Subroutine Funcs(Np,Gp,x,y,F,dFx,dFy,G,dGx,dGy,E,iFlag)
      Dimension E(728)
      iFlag=0
*Note:  While compiling in DOS platform overflow checking commands
*   may need to be changed
      F=(-1.0)*Float(Np)
      G=(-1.0)*2.0/Gp
      dFx=0.0
      dFy=0.0
      dGx=0.0
      dGy=0.0
      Do 10 I=1,Np
        Dr =((E(I)-y)**2+x**2)**0.5
        dDr=Dr*((E(I)-y)**2+x**2)
        F=F+(1.0 - (E(I)-y)/Dr)
        dFx=dFx + x*(E(I)-y)/dDr
        dFy=dFy + x**2/dDr
        G=G+1.0/Dr
        dGx=dGx - x/dDr
        dGy=dFy + (E(I)-y)/dDr
10    Continue
      iFlag=1
      End

