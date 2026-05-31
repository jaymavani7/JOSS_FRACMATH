!=======================================================================
!  CDM_UMAT_2D.FOR
!  Abaqus/Standard UMAT for 2D plane-stress isotropic scalar continuum
!  damage mechanics. Matches the MATLAB material law; Abaqus uses CELENT
!
!  Constitutive law:
!     - Modified von Mises equivalent strain (de Vree et al. 1995)
!     - Exponential softening
!     - Crack-band regularization using Abaqus's CELENT
!
!  Element type: CPS3 (3-node linear plane-stress triangle, 1 IP).
!  Also works with CPS4 / CPS4R.
!  Tangent returned is the SECANT tangent (1 - omega)*C0.
!
!  Material constants (PROPS), 5 values:
!     PROPS(1) = E        Young's modulus       [MPa]
!     PROPS(2) = nu       Poisson ratio
!     PROPS(3) = ft       tensile strength      [MPa]
!     PROPS(4) = GF       mode-I fracture energy [N/mm]
!     PROPS(5) = fc_ft    fc / ft  (10 for concrete)
!
!  State variables (STATEV), 2 values:
!     STATEV(1) = kappa   history variable
!     STATEV(2) = omega   damage 0..1
!
!  Default material (Gregoire concrete, matches MATLAB):
!     E = 37000 MPa,  nu = 0.20
!     ft = 3.5 MPa,   GF = 0.090 N/mm
!     fc/ft = 10
!
!  Authors: [Name to be added]
!=======================================================================
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,
     1  RPL, DDSDDT, DRPLDE, DRPLDT,
     2  STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,
     3  NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,
     4  CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, JSTEP, KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS),
     1  DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2  TIME(2), PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3),
     3  DROT(3,3), DFGRD0(3,3), DFGRD1(3,3), JSTEP(4)

      DIMENSION EPS(3), C0(3,3)
      REAL*8 E, ANU, FT, GF, FCFT, EPS0
      REAL*8 KAPPA, OMEGA, OMEGA_NEW
      REAL*8 EXX, EYY, GXY
      REAL*8 EM, RAD, E1, E2, E3
      REAL*8 I1, J2, A1, A2, A3, A4, INSIDE, EQ
      REAL*8 EF, H, BUF, FAC
      REAL*8 OMAX
      PARAMETER (OMAX = 0.999999999999D0)

!-------- 1. Read material constants ----------------------------------
      E    = PROPS(1)
      ANU  = PROPS(2)
      FT   = PROPS(3)
      GF   = PROPS(4)
      FCFT = PROPS(5)
      EPS0 = FT / E

!-------- 2. Read history --------------------------------------------
      KAPPA = STATEV(1)
      OMEGA = STATEV(2)

!-------- 3. Total strain (plane stress, NTENS = 3) ------------------
!     Abaqus convention: STRAN = [exx, eyy, gxy], engineering shear.
      EPS(1) = STRAN(1) + DSTRAN(1)
      EPS(2) = STRAN(2) + DSTRAN(2)
      EPS(3) = STRAN(3) + DSTRAN(3)
      EXX = EPS(1)
      EYY = EPS(2)
      GXY = EPS(3)

!-------- 4. Principal in-plane strains + plane-stress ezz -----------
      EM  = 0.5D0 * (EXX + EYY)
      RAD = SQRT(0.25D0*(EXX - EYY)**2 + 0.25D0*GXY*GXY)
      E1  = EM + RAD
      E2  = EM - RAD
      E3  = -(ANU / (1.0D0 - ANU)) * (E1 + E2)

      I1 = E1 + E2 + E3
      J2 = 0.5D0 * ((E1-E2)**2 + (E2-E3)**2 + (E3-E1)**2)

!-------- 5. Modified von Mises equivalent strain --------------------
      A1 = (FCFT - 1.0D0) / (2.0D0 * FCFT * (1.0D0 - 2.0D0*ANU))
      A2 = 1.0D0 / (2.0D0 * FCFT)
      A3 = ((FCFT - 1.0D0) / (1.0D0 - 2.0D0*ANU))**2
      A4 = 12.0D0 * FCFT / (1.0D0 + ANU)**2

      INSIDE = A3 * I1*I1 + A4 * J2
      IF (INSIDE .LT. 0.0D0) INSIDE = 0.0D0
      EQ = A1*I1 + A2*SQRT(INSIDE)
      IF (EQ .LT. 0.0D0) EQ = 0.0D0

!-------- 6. Update history (Kuhn-Tucker) ----------------------------
      IF (EQ .GT. KAPPA) KAPPA = EQ

!-------- 7. Exponential softening with crack-band scaling -----------
!     ef = eps0/2 + GF / (h * ft),  h = CELENT (Abaqus characteristic
!     element length).
      H = CELENT
      BUF = GF / (H * FT)
      EF = 0.5D0*EPS0 + BUF
      IF (EF .LE. EPS0) EF = EPS0 + 1.0D-12

      IF (KAPPA .LE. EPS0) THEN
         OMEGA_NEW = 0.0D0
      ELSE
         OMEGA_NEW = 1.0D0 - (EPS0/KAPPA) *
     &       EXP(-(KAPPA - EPS0) / (EF - EPS0))
      END IF

      IF (OMEGA_NEW .LT. OMEGA) OMEGA_NEW = OMEGA   ! monotonic
      IF (OMEGA_NEW .LT. 0.0D0) OMEGA_NEW = 0.0D0
      IF (OMEGA_NEW .GT. OMAX)  OMEGA_NEW = OMAX
      OMEGA = OMEGA_NEW

!-------- 8. Plane-stress elasticity matrix --------------------------
      DO I = 1, NTENS
         DO J = 1, NTENS
            C0(I,J) = 0.0D0
         END DO
      END DO
      C0(1,1) = E / (1.0D0 - ANU*ANU)
      C0(1,2) = E * ANU / (1.0D0 - ANU*ANU)
      C0(2,1) = C0(1,2)
      C0(2,2) = C0(1,1)
      C0(3,3) = E / (2.0D0 * (1.0D0 + ANU))

!-------- 9. Stress = (1 - omega) * C0 * eps -------------------------
      FAC = 1.0D0 - OMEGA
      DO I = 1, NTENS
         STRESS(I) = 0.0D0
         DO J = 1, NTENS
            STRESS(I) = STRESS(I) + FAC * C0(I,J) * EPS(J)
         END DO
         DO J = 1, NTENS
            DDSDDE(I,J) = FAC * C0(I,J)
         END DO
      END DO

!-------- 10. Store updated state -----------------------------------
      STATEV(1) = KAPPA
      STATEV(2) = OMEGA

      RETURN
      END
