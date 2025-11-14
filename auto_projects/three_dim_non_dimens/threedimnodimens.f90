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

  DOUBLE PRECISION v, y, z, I, gs, Vs, Vna, gk, Vk, gl, Vl, eps, a, m, alpha_y, beta_y, beta_z, alpha_z, kv

  DOUBLE PRECISION MINF, YINF, ZINF, TAUY, TAUZ

  MINF(v) = 0.5*(1+DTANH((kv*v-m)/a))
  YINF(v) = 0.5*(1+DTANH((kv*v - beta_y)/alpha_y))
  ZINF(v) = 0.5*(1+DTANH((kv*v - beta_z)/alpha_z))
  TAUY(v) = (DCOSH((kv*v - beta_y)/(2*alpha_y)))
  TAUZ(v) = (DCOSH((kv*v - beta_z)/(2*alpha_z)))

  v = U(1)
  y = U(2)
  z = U(3)

  I = PAR(1)
  gs = PAR(2)
  Vs = PAR(3)
  Vna = PAR(4)
  gk = PAR(5)
  Vk = PAR(6)
  gl = PAR(7)
  Vl = PAR(8)
  eps = PAR(9)
  a = PAR(10)
  m = PAR(16)
  alpha_y = PAR(17)
  beta_y = PAR(18)
  beta_z = PAR(19)
  alpha_z = PAR(20)
  kv = PAR(24)

  F(1) = I-MINF(v)*(v-Vna)-gs*z*(v-Vs)-gk*y*(v-Vk)-gl*(v-Vl)
  F(2) = eps*(YINF(v)-y)*TAUY(v)
  F(3) = eps*(ZINF(v)-z)*TAUZ(v)

END SUBROUTINE FUNC

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE STPNT(NDIM,U,PAR,T)
!--------- -----

! Input arguments :
!      NDIM   :   Dimension of the algebraic or ODE system 

! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values

! Note : For time- or space-dependent solutions this subroutine has
!        the scalar input parameter T contains the varying time or space
!        variable value.
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

! Initialize the equation parameters
  PAR(1) = 0
  PAR(2) = 0.15d0
  PAR(3) = 0.5d0
  PAR(4) = 0.5d0
  PAR(5) = 1
  PAR(6) = -1
  PAR(7) = 0.1d0
  PAR(8) = -0.7d0
  PAR(9) = 0.0015d0
  PAR(10) = 18
  PAR(16) = -1.2d0
  PAR(17) = 10
  PAR(18) = -10
  PAR(19) = -21
  PAR(20) = 15
  PAR(24) = 100

! Initialize the solution
  U(1) = -0.6907786130905151d0
  U(2) = 7.388543963315897e-06
  U(3) = 0.001641697133891284d0
   
END SUBROUTINE STPNT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!--------- ----

! Boundary Conditions

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      PAR    :   Equation parameters
!      ICP    :   Array indicating the free parameter(s)
!      NBC    :   Number of boundary conditions
!      U0     :   State variable values at the left boundary
!      U1     :   State variable values at the right boundary

! Values to be returned :
!      FB     :   The values of the boundary condition functions 

! Normally unused Jacobian arguments : IJAC, DBC (see manual)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
  DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

!X FB(1)=
!X FB(2)=

END SUBROUTINE BCND

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
!--------- ----

! Integral Conditions

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      PAR    :   Equation parameters
!      ICP    :   Array indicating the free parameter(s)
!      NINT   :   Number of integral conditions
!      U      :   Value of the vector function U at `time' t

! The following input arguments, which are normally not needed,
! correspond to the preceding point on the solution branch
!      UOLD   :   The state vector at 'time' t
!      UDOT   :   Derivative of UOLD with respect to arclength
!      UPOLD  :   Derivative of UOLD with respect to `time'

! Normally unused Jacobian arguments : IJAC, DINT

! Values to be returned :
!      FI     :   The value of the vector integrand 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
  DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

!X FI(1)=

END SUBROUTINE ICND

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE FOPT(NDIM,U,ICP,PAR,IJAC,FS,DFDU,DFDP)
!--------- ----
!
! Defines the objective function for algebraic optimization problems
!
! Supplied variables :
!      NDIM   :   Dimension of the state equation
!      U      :   The state vector
!      ICP    :   Indices of the control parameters
!      PAR    :   The vector of control parameters
!
! Values to be returned :
!      FS      :   The value of the objective function
!
! Normally unused Jacobian argument : IJAC, DFDP

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: FS
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM),DFDP(*)

!X FS=

END SUBROUTINE FOPT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE PVLS(NDIM,U,PAR)
!--------- ----

  IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      INTEGER NDX,NCOL,NTST
!---------------------------------------------------------------------- 
! NOTE : 
! Parameters set in this subroutine should be considered as ``solution 
! measures'' and be used for output purposes only.
! 
! They should never be used as `true'' continuation parameters. 
!
! They may, however, be added as ``over-specified parameters'' in the 
! parameter list associated with the AUTO-Constant NICP, in order to 
! print their values on the screen and in the ``p.xxx file.
!
! They may also appear in the list associated with AUTO-constant NUZR.
!
!---------------------------------------------------------------------- 
! For algebraic problems the argument U is, as usual, the state vector.
! For differential equations the argument U represents the approximate 
! solution on the entire interval [0,1]. In this case its values can
! be accessed indirectly by calls to GETP, as illustrated below, or
! by obtaining NDIM, NCOL, NTST via GETP and then dimensioning U as
! U(NDIM,0:NCOL*NTST) in a seperate subroutine that is called by PVLS.
!---------------------------------------------------------------------- 

! Set PAR(12) equal to the minimum of U(1)
       PAR(12)=GETP('MIN',1,U)

! Set PAR(13) equal to the stability at the point
       PAR(13)=GETP('STA', 1, U)

END SUBROUTINE PVLS
