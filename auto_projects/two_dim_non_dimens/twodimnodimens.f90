!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   twodim.f90 - adapted from Presscot XPP .ode files
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

  DOUBLE PRECISION v, w, I, beta, eps, Vna, gs, Vs, gl, Vl, a, m, alpha
  DOUBLE PRECISION MINF, WINF, TAUW
  
  MINF(V) = 0.5*(1+DTANH(a*v-m))
  WINF(V) = 0.5*(1+DTANH(alpha*v-beta))
  TAUW(V) = (DCOSH(0.5*(alpha*v-beta)))

  v = U(1)
  w = U(2)

  I = PAR(1)
  beta = PAR(2)
  eps = PAR(3)
  Vna = PAR(4)
  gs = PAR(5)
  Vs = PAR(6)
  gl = PAR(7)
  Vl = PAR(8)
  a = PAR(9)
  m = PAR(16)
  alpha = PAR(17)

  F(1) = I-MINF(v)*(v-Vna)-gs*w*(v-Vs)-gl*(v-Vl)
  F(2) = eps*(WINF(v)-w)*TAUW(v)

END SUBROUTINE FUNC

SUBROUTINE STPNT(NDIM,U,PAR,T)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

! Initialize the equation parameters
  PAR(1) = 0.0              ! I 
  PAR(2) = 0.0              ! beta = beta_w / gamma_w
  PAR(3) = REAL(0.15*2)/20  ! eps = phi_w *c / g_fast
  PAR(4) = REAL(50)/100     ! Vna = Vna / kv
  PAR(5) = REAL(20)/20      ! gs = g_slow / g_fast 
  PAR(6) = REAL(-100)/100   ! Vs = V_slow / kv
  PAR(7) = REAL(2)/20       ! gl = g_leak / g_fast
  PAR(8) = REAL(-70)/100    ! Vl = V_leak / kv
  PAR(9) = REAL(100)/18     ! a = kv / gamma_m
  PAR(16) = REAL(-1.2)/18   ! m = beta_m / gamma_m 
  PAR(17) = REAL(100)/10    ! alpha = kv / gamma_w

! Initialize the solution
  U(1) = -0.6938893890380859
  U(2) = 9.396213158652245e-07
  
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

! Set PAR(12) equal to the minimum of U(1)
       PAR(12)=GETP('MIN',1,U)

! Set PAR(13) equal to the stability at the point
       PAR(13)=GETP('STA', 1, U)

END SUBROUTINE PVLS
