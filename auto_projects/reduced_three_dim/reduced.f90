!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   twodim.f90 - adapted from Presscot XPP .ode files
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!--------- ----

! Evaluates the algebraic equations or ODE right hand side

! Input arguments :
!      NDIM   :   Dimension of the algebraic or ODE system 
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters

! Values to be returned :
!      F      :   Equation or ODE right hand side values

! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

  DOUBLE PRECISION V, z, i_dc, gsub, Vsub, Vna, gk, Vk, gl, Vl, C, phi_y, beta_m, gamma_m, gamma_y, beta_y, beta_z, gamma_z, gna, phi_z

  DOUBLE PRECISION minf, yinf, zinf, phi, tauy, tauz

  minf(V)=.5*(1+tanh((V-beta_m)/gamma_m))

  phi(V, z) = (i_dc-gna*minf(V)*(V-Vna)-gl*(V-Vl)-gsub*z*(V-Vsub))/(gk*(V-VK))
  yinf(V)=.5*(1+TANH((V-beta_y)/gamma_y))
  tauy(V)=1/COSH((V-beta_y)/(2*gamma_y))

  zinf(V)=.5*(1+TANH((V-beta_z)/gamma_z))
  tauz(V)=1/COSH((V-beta_z)/(2*gamma_z))

  V = U(1)
  z = U(2)

  i_dc = PAR(1)       !0.0d0
  gsub = PAR(2)       !8.0d0
  Vsub = PAR(3)       !-100.0d0
  Vna = PAR(4)        !50.0d0
  gk = PAR(5)         !20.0d0
  Vk = PAR(6)         !-100.0d0
  gl = PAR(7)         !2.0d0
  Vl = PAR(8)         !-70.0d0
  c = PAR(9)          !2.0d0
  phi_y = PAR(10)     !0.0015d0
  beta_m = PAR(16)    !-1.2d0
  gamma_m = PAR(17)   !18.0d0
  gamma_y = PAR(18)   !10.0d0
  beta_y = PAR(19)    !-10.0d0
  beta_z = PAR(20)    !-21.0d0
  gamma_z = PAR(24)   !15.0d0
  gna = PAR(25)       !20.0d0
  phi_z = PAR(26)     !0.0015d0

  F(1) = ((-gsub*(V-Vsub))/c)*(phi_z*(zinf(V)-z)/tauz(V)) + ((-gk*(V-VK))/c)*(phi_y*(yinf(V)-phi(V, z))/tauy(V))
  F(2) = - (((-gna*minf(V) -gk*phi(V,z) -gl -gsub*z -gna*(V-Vna)*(.5*(1/(cosh((V-beta_m)/gamma_m))**2)/gamma_m)))/c) * (phi_z*(zinf(V)-z)/tauz(V))

END SUBROUTINE FUNC

!-----------------------------------------------------------------------

SUBROUTINE STPNT(NDIM,U,PAR,T)
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

! Initialize the equation parameters
  PAR(1:10) = (/0.0d0, 8.0d0, -100.0d0, 50.0d0, 20.0d0, -100.0d0, 2.0d0, -70.0d0, 2.0d0, 0.0015d0 /)
  PAR(16:20) = (/ -1.2d0, 18.0d0, 10.0d0,-10.0d0, -21.0d0 /)
  PAR(24) = 15.0d0
  PAR(25) = 20.0d0 
  PAR(26) = 0.0015d0

! Initialize the solution
  U(1) = -69.59233093261719d0
  U(2) = 0.001554027083329856d0
   
END SUBROUTINE STPNT

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
END SUBROUTINE BCND

SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
END SUBROUTINE ICND

SUBROUTINE FOPT(NDIM,U,ICP,PAR,IJAC,FS,DFDU,DFDP)
END SUBROUTINE FOPT

SUBROUTINE PVLS(NDIM,U,PAR)


  IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      INTEGER NDX,NCOL,NTST

! Set PAR(13) equal to the stability at the point
       PAR(12)=GETP('STA', 1, U)

END SUBROUTINE PVLS