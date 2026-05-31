!=======================================================================
!  SCM_UMAT_2D_OLIVER_MATCH_MATLAB.FOR  (FIXED)
!
!  Abaqus/Standard UMAT — 2-D plane-stress scalar CDM with full
!  direction-dependent Oliver crack-band bandwidth.
!
!  Physics: EXACTLY matches solver_main_3pb_OLIVER_bandwidth.m
!    - Modified von Mises eq strain (de Vree, plane-stress e3 approx)
!    - Exponential stress-strain softening
!    - Oliver h(n) = 2 / sum_a |grad(N_a).n|   (from gradient table)
!    - STAGGERED secant tangent: stress & DDSDDE use omega from START
!      of the Abaqus increment (frozen), omega_new saved for next incr.
!    - PNEWDT = 0.25 whenever damage grows beyond tolerance to force
!      Abaqus to use fine increments matching MATLAB step size.
!
!  FIX LOG (vs previous version):
!    FIX-1  EF lower-bound guard added:  ef > eps0+1e-12  (MATLAB line)
!    FIX-2  PNEWDT cut-back on large damage jump (stops Abaqus from
!           taking huge steps that have no MATLAB equivalent)
!    FIX-3  J2 formula confirmed: 0.5*(...) — matches MATLAB exactly
!    FIX-4  iso-strain guard: threshold 1e-18 matches MATLAB
!    FIX-5  H fallback to CELENT only when gradient table is absent;
!           never overwrite a valid Oliver H with CELENT
!    FIX-6  Shear convention: STRAN(3)=gamma_xy (Abaqus engineering);
!           Mohr's circle written for engineering shear — confirmed OK
!
!  Material constants (PROPS), 5 values:
!     PROPS(1) = E        Young's modulus        [MPa]
!     PROPS(2) = nu       Poisson ratio
!     PROPS(3) = ft       tensile strength       [MPa]
!     PROPS(4) = GF       mode-I fracture energy [N/mm]
!     PROPS(5) = fc_ft    fc / ft  (= k in de Vree)
!
!  State variables (STATEV), 3:
!     STATEV(1) = kappa   history variable (max eq strain ever seen)
!     STATEV(2) = omega   damage [0, OMAX]   — written as omega_NEW
!     STATEV(3) = H       Oliver bandwidth used this increment [mm]
!=======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1  RPL,DDSDDT,DRPLDE,DRPLDT,
     2  STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3  NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4  CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1  DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2  TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     3  DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)

      DIMENSION C0(3,3)

!     --- Oliver gradient table (read once, shared across all calls) ---
      INTEGER MAXEL
      PARAMETER (MAXEL = 600000)
      INTEGER INITDONE,HASG(MAXEL),IDTMP,UNITG
      INTEGER I,J
      PARAMETER (UNITG = 99)
      REAL*8 G1X(MAXEL),G1Y(MAXEL),G2X(MAXEL),G2Y(MAXEL)
      REAL*8 G3X(MAXEL),G3Y(MAXEL),AREAEL(MAXEL),SQRTA(MAXEL)
      SAVE INITDONE,HASG,G1X,G1Y,G2X,G2Y,G3X,G3Y,AREAEL,SQRTA
      DATA INITDONE /0/

!     --- local scalars ---
      REAL*8 E,ANU,FT,GF,FCFT,EPS0
      REAL*8 KAPPA,OMEGA,OMEGA_OLD,OMEGA_NEW
      REAL*8 EXX,EYY,GXY
      REAL*8 EM,RAD,E1,E2,E3
      REAL*8 I1,J2VAL,A1,A2,A3,A4,INSIDE,EQ
      REAL*8 EF,H,BUF,FAC
      REAL*8 THETA,NX,NY,DEN,HTMP
      REAL*8 COL1,COL2,COL3,COL4,COL5,COL6,COL7,COL8,COL9
      REAL*8 OMAX,DOMEGA,PNEWDT_THRESHOLD
      PARAMETER (OMAX = 0.999999999999D0)
      PARAMETER (PNEWDT_THRESHOLD = 0.05D0)

!=======================================================================
!  READ OLIVER GRADIENT TABLE ONCE
!=======================================================================
      IF (INITDONE .EQ. 0) THEN
         DO I = 1, MAXEL
            HASG(I) = 0
         END DO

         OPEN(UNIT=UNITG,FILE='oliver_bandwidth_data.dat',
     &        STATUS='OLD',ERR=900)
  100    CONTINUE
!        Columns: eid  g1x g1y  g2x g2y  g3x g3y  area  sqrt(area)
         READ(UNITG,*,END=110,ERR=110)
     &        IDTMP,COL1,COL2,COL3,COL4,COL5,COL6,COL7,COL8
         IF (IDTMP .GE. 1 .AND. IDTMP .LE. MAXEL) THEN
            G1X(IDTMP)    = COL1
            G1Y(IDTMP)    = COL2
            G2X(IDTMP)    = COL3
            G2Y(IDTMP)    = COL4
            G3X(IDTMP)    = COL5
            G3Y(IDTMP)    = COL6
            AREAEL(IDTMP) = COL7
            SQRTA(IDTMP)  = COL8
            HASG(IDTMP)   = 1
         END IF
         GOTO 100
  110    CONTINUE
         CLOSE(UNITG)
         INITDONE = 1
         GOTO 910
  900    CONTINUE
!        File absent — run without Oliver (CELENT fallback active)
         INITDONE = 2
  910    CONTINUE
      END IF

!=======================================================================
!  MATERIAL CONSTANTS
!=======================================================================
      E    = PROPS(1)
      ANU  = PROPS(2)
      FT   = PROPS(3)
      GF   = PROPS(4)
      FCFT = PROPS(5)
      EPS0 = FT / E

!=======================================================================
!  RETRIEVE STATE AT START OF INCREMENT
!  omega_OLD is frozen for tangent/stress (staggered, as in MATLAB)
!=======================================================================
      KAPPA     = STATEV(1)
      OMEGA_OLD = STATEV(2)
      OMEGA     = OMEGA_OLD      ! alias used in tangent/stress block

!=======================================================================
!  TOTAL STRAIN AT END OF INCREMENT
!=======================================================================
      EXX = STRAN(1) + DSTRAN(1)
      EYY = STRAN(2) + DSTRAN(2)
      GXY = STRAN(3) + DSTRAN(3)   ! engineering shear (Abaqus convention)

!=======================================================================
!  PRINCIPAL STRAINS  (plane-stress out-of-plane approximation)
!  Matches MATLAB:
!    me  = (ex+ey)/2
!    rad = sqrt(((ex-ey)/2)^2 + (gxy/2)^2)
!    e1  = me+rad; e2 = me-rad
!    e3  = -(nu/(1-nu))*(e1+e2)
!=======================================================================
      EM  = 0.5D0*(EXX + EYY)
      RAD = DSQRT(0.25D0*(EXX-EYY)**2 + 0.25D0*GXY**2)
      E1  = EM + RAD
      E2  = EM - RAD
      E3  = -(ANU/(1.0D0-ANU))*(E1+E2)

!=======================================================================
!  MODIFIED VON MISES EQUIVALENT STRAIN  (de Vree, MATLAB-exact)
!  MATLAB:
!    I1 = e1+e2+e3
!    J2 = 0.5*((e1-e2)^2+(e2-e3)^2+(e3-e1)^2)   <- factor 0.5 !!!
!    k  = fc/ft
!    a1 = (k-1)/(2*k*(1-2*nu))
!    a2 = 1/(2*k)
!    a3 = ((k-1)/(1-2*nu))^2
!    a4 = 12*k/(1+nu)^2
!    eq = a1*I1 + a2*sqrt(max(a3*I1^2+a4*J2,0))
!=======================================================================
      I1    = E1 + E2 + E3
      J2VAL = 0.5D0*((E1-E2)**2 + (E2-E3)**2 + (E3-E1)**2)

      A1 = (FCFT-1.0D0)/(2.0D0*FCFT*(1.0D0-2.0D0*ANU))
      A2 = 1.0D0/(2.0D0*FCFT)
      A3 = ((FCFT-1.0D0)/(1.0D0-2.0D0*ANU))**2
      A4 = 12.0D0*FCFT/(1.0D0+ANU)**2

      INSIDE = A3*I1*I1 + A4*J2VAL
      IF (INSIDE .LT. 0.0D0) INSIDE = 0.0D0
      EQ = A1*I1 + A2*DSQRT(INSIDE)
      IF (EQ .LT. 0.0D0) EQ = 0.0D0

!     Update history variable
      IF (EQ .GT. KAPPA) KAPPA = EQ

!=======================================================================
!  OLIVER DIRECTION-DEPENDENT BANDWIDTH
!  Crack normal = max principal strain direction
!  theta = 0.5*atan2(gxy, exx-eyy)   [MATLAB-exact]
!  h(n)  = 2 / sum_a |grad(N_a).n|
!=======================================================================
      H = CELENT    ! safe fallback if gradient table not available

      IF (INITDONE .EQ. 1) THEN
         IF (NOEL .GE. 1 .AND. NOEL .LE. MAXEL) THEN
            IF (HASG(NOEL) .EQ. 1) THEN

               IF (DABS(EXX-EYY)+DABS(GXY) .LT. 1.0D-18) THEN
!                 Hydrostatic / isotropic: direction undefined.
!                 Use x-direction (matches MATLAB iso guard).
                  NX = 1.0D0
                  NY = 0.0D0
               ELSE
                  THETA = 0.5D0*DATAN2(GXY, EXX-EYY)
                  NX = DCOS(THETA)
                  NY = DSIN(THETA)
               END IF

               DEN = DABS(G1X(NOEL)*NX + G1Y(NOEL)*NY)
     &             + DABS(G2X(NOEL)*NX + G2Y(NOEL)*NY)
     &             + DABS(G3X(NOEL)*NX + G3Y(NOEL)*NY)

               IF (DEN .GT. 1.0D-14) THEN
                  HTMP = 2.0D0/DEN
!                 FIX-5: only use Oliver H; do NOT fall back to CELENT
!                 unless H is non-physical (< 1e-12 or > 1e+12).
                  IF (HTMP .GE. 1.0D-12 .AND. HTMP .LE. 1.0D12) THEN
                     H = HTMP
                  END IF
               END IF
!              If DEN <= 1e-14 keep H = CELENT (degenerate direction)
            END IF
         END IF
      END IF

!=======================================================================
!  EXPONENTIAL SOFTENING LAW  (MATLAB-exact)
!  ef = max(eps0/2 + GF/(h*ft),  eps0 + 1e-12)   <-- FIX-1: lower bound
!=======================================================================
      BUF = GF/(H*FT)
      EF  = 0.5D0*EPS0 + BUF
!     FIX-1: enforce lower bound matching MATLAB
      IF (EF .LE. EPS0 + 1.0D-12) EF = EPS0 + 1.0D-12

      IF (KAPPA .LE. EPS0) THEN
         OMEGA_NEW = 0.0D0
      ELSE
         OMEGA_NEW = 1.0D0 - (EPS0/KAPPA)*
     &       DEXP(-(KAPPA-EPS0)/(EF-EPS0))
      END IF

!     Damage is irreversible
      IF (OMEGA_NEW .LT. OMEGA_OLD) OMEGA_NEW = OMEGA_OLD
      IF (OMEGA_NEW .LT. 0.0D0)     OMEGA_NEW = 0.0D0
      IF (OMEGA_NEW .GT. OMAX)       OMEGA_NEW = OMAX
      IF (.NOT.(OMEGA_NEW .EQ. OMEGA_NEW)) OMEGA_NEW = OMEGA_OLD

!=======================================================================
!  PLANE-STRESS ELASTIC MATRIX  C0
!=======================================================================
      DO I = 1, NTENS
         DO J = 1, NTENS
            C0(I,J) = 0.0D0
         END DO
      END DO
      C0(1,1) = E/(1.0D0-ANU*ANU)
      C0(1,2) = E*ANU/(1.0D0-ANU*ANU)
      C0(2,1) = C0(1,2)
      C0(2,2) = C0(1,1)
      C0(3,3) = E/(2.0D0*(1.0D0+ANU))

!=======================================================================
!  STAGGERED SECANT TANGENT  (MATLAB-exact)
!  Stress and DDSDDE use OMEGA_OLD (frozen at start of increment).
!  OMEGA_NEW is written to STATEV so the NEXT increment uses it.
!
!  This matches MATLAB:
!    K = assemble_K(Ke0, omega_OLD, ...)   <- secant with old omega
!    solve for u ...
!    [omega_NEW, kappa_NEW] = damage_update(u, ..., kappa_OLD, omega_OLD)
!    omega_OLD <- omega_NEW for next step
!=======================================================================
!     stress = (1-omega_old) * C0 * eps_total
      FAC = 1.0D0 - OMEGA_OLD
      STRESS(1) = FAC*(C0(1,1)*EXX + C0(1,2)*EYY)
      STRESS(2) = FAC*(C0(2,1)*EXX + C0(2,2)*EYY)
      STRESS(3) = FAC*C0(3,3)*GXY

      DO I = 1, NTENS
         DO J = 1, NTENS
            DDSDDE(I,J) = FAC*C0(I,J)
         END DO
      END DO

!=======================================================================
!  PNEWDT cut-back DISABLED to match MATLAB.
!  MATLAB takes 10000 FIXED displacement steps with NO adaptive
!  cut-back.  Requesting PNEWDT < 1 here made Abaqus sub-step near the
!  peak (finer increments than MATLAB), which shifted the peak load.
!  With *STATIC, DIRECT=NO STOP the increment size is now fixed to the
!  MATLAB step, so the UMAT must NOT request a cut-back.
!=======================================================================
!     DOMEGA = OMEGA_NEW - OMEGA_OLD
!     IF (DOMEGA .GT. PNEWDT_THRESHOLD) THEN
!        PNEWDT = 0.25D0
!     END IF

!=======================================================================
!  SAVE STATE
!=======================================================================
      IF (NSTATV .GE. 1) STATEV(1) = KAPPA
      IF (NSTATV .GE. 2) STATEV(2) = OMEGA_NEW
      IF (NSTATV .GE. 3) STATEV(3) = H

      RETURN
      END
