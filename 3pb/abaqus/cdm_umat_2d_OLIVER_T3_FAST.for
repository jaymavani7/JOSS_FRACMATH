C=======================================================================
C  CDM_UMAT_2D_OLIVER_T3_FAST.FOR
C
C  Faster Abaqus/Standard UMAT for CPS3 plane-stress scalar CDM.
C  Uses MATLAB-equivalent Oliver direction-dependent T3 crack-band width:
C
C       h(n) = 2 / SUM_a | grad(N_a) dot n |,  a = 1..3
C
C  Required job-folder file:
C       oliver_t3_gradN.dat
C  Format:
C       NOEL, g1x, g1y, g2x, g2y, g3x, g3y
C
C  Speed notes:
C    - UMAT is still called point-by-point by Abaqus/Standard; true
C      block vectorization is not possible inside a Standard UMAT.
C    - This file removes small matrix loops in the stress/tangent update.
C    - Use *DEPVAR, n=2 for fastest output: STATEV(1)=kappa, STATEV(2)=omega.
C      If n>=4, STATEV(3)=h and STATEV(4)=Oliver flag are stored for checking.
C=======================================================================
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,
     1  RPL, DDSDDT, DRPLDE, DRPLDT,
     2  STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED,
     3  CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS,
     4  DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     5  KSPT, JSTEP, KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      REAL*8 SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, DTEMP
      REAL*8 PNEWDT, CELENT
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS),
     1  DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2  TIME(2), PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3),
     3  DROT(3,3), DFGRD0(3,3), DFGRD1(3,3), JSTEP(4)

      REAL*8 E, ANU, FT, GF, FCFT, EPS0
      REAL*8 KAPPA, OMEGA, OMEGA_NEW, FAC
      REAL*8 EXX, EYY, GXY, EM, RAD, E1, E2, E3
      REAL*8 I1, J2, A1, A2, A3, A4, INSIDE, EQ
      REAL*8 EF, H, C11, C12, C33, DENOM
      INTEGER I, J, IHFLAG
      REAL*8 OMAX
      PARAMETER (OMAX = 0.999999999999D0)

C----- Material constants ---------------------------------------------
      E    = PROPS(1)
      ANU  = PROPS(2)
      FT   = PROPS(3)
      GF   = PROPS(4)
      FCFT = PROPS(5)
      EPS0 = FT / E

C----- State variables -------------------------------------------------
      KAPPA = STATEV(1)
      OMEGA = STATEV(2)

C----- Total strain; Abaqus uses engineering shear gamma_xy -----------
      EXX = STRAN(1) + DSTRAN(1)
      EYY = STRAN(2) + DSTRAN(2)
      GXY = STRAN(3) + DSTRAN(3)

C----- Principal strains and plane-stress out-of-plane strain ---------
      EM  = 0.5D0 * (EXX + EYY)
      RAD = DSQRT(0.25D0*(EXX - EYY)*(EXX - EYY)
     &          + 0.25D0*GXY*GXY)
      E1  = EM + RAD
      E2  = EM - RAD
      E3  = -(ANU / (1.0D0 - ANU)) * (E1 + E2)

      I1 = E1 + E2 + E3
      J2 = 0.5D0 * ((E1-E2)*(E1-E2) + (E2-E3)*(E2-E3)
     &             + (E3-E1)*(E3-E1))

C----- Modified von Mises equivalent strain ---------------------------
      A1 = (FCFT - 1.0D0) / (2.0D0 * FCFT *
     &     (1.0D0 - 2.0D0*ANU))
      A2 = 1.0D0 / (2.0D0 * FCFT)
      A3 = ((FCFT - 1.0D0) / (1.0D0 - 2.0D0*ANU))**2
      A4 = 12.0D0 * FCFT / ((1.0D0 + ANU)*(1.0D0 + ANU))

      INSIDE = A3 * I1*I1 + A4 * J2
      IF (INSIDE .LT. 0.0D0) INSIDE = 0.0D0
      EQ = A1*I1 + A2*DSQRT(INSIDE)
      IF (EQ .LT. 0.0D0) EQ = 0.0D0
      IF (EQ .GT. KAPPA) KAPPA = EQ

C----- Oliver bandwidth from table; CELENT only if table unavailable ---
      CALL OLIVER_H_T3_FAST(NOEL, EXX, EYY, GXY, CELENT, H, IHFLAG)

      IF (NSTATV .GE. 3) STATEV(3) = H
      IF (NSTATV .GE. 4) STATEV(4) = DBLE(IHFLAG)

C----- Exponential softening, same as MATLAB --------------------------
      EF = 0.5D0*EPS0 + GF/(H*FT)
      IF (EF .LE. EPS0) EF = EPS0 + 1.0D-12

      IF (KAPPA .LE. EPS0) THEN
         OMEGA_NEW = 0.0D0
      ELSE
         OMEGA_NEW = 1.0D0 - (EPS0/KAPPA) *
     &       DEXP(-(KAPPA - EPS0) / (EF - EPS0))
      END IF

      IF (OMEGA_NEW .LT. OMEGA) OMEGA_NEW = OMEGA
      IF (OMEGA_NEW .LT. 0.0D0) OMEGA_NEW = 0.0D0
      IF (OMEGA_NEW .GT. OMAX)  OMEGA_NEW = OMAX
      OMEGA = OMEGA_NEW

C----- Unrolled plane-stress stress and secant tangent ----------------
      DENOM = 1.0D0 - ANU*ANU
      C11 = E / DENOM
      C12 = E * ANU / DENOM
      C33 = E / (2.0D0 * (1.0D0 + ANU))
      FAC = 1.0D0 - OMEGA

      STRESS(1) = FAC * (C11*EXX + C12*EYY)
      STRESS(2) = FAC * (C12*EXX + C11*EYY)
      STRESS(3) = FAC * (C33*GXY)

      DO I = 1, NTENS
         DO J = 1, NTENS
            DDSDDE(I,J) = 0.0D0
         END DO
      END DO
      DDSDDE(1,1) = FAC*C11
      DDSDDE(1,2) = FAC*C12
      DDSDDE(2,1) = FAC*C12
      DDSDDE(2,2) = FAC*C11
      DDSDDE(3,3) = FAC*C33

C----- Store state -----------------------------------------------------
      STATEV(1) = KAPPA
      STATEV(2) = OMEGA

      RETURN
      END

C=======================================================================
C  Oliver bandwidth for CPS3 element NOEL.
C=======================================================================
      SUBROUTINE OLIVER_H_T3_FAST(NOEL, EXX, EYY, GXY, CELENT,
     &                            H, IFLAG)

      INCLUDE 'ABA_PARAM.INC'

      INTEGER NOEL, IFLAG
      LOGICAL HAVE
      REAL*8 EXX, EYY, GXY, CELENT, H
      REAL*8 G1X, G1Y, G2X, G2Y, G3X, G3Y
      REAL*8 THETA, NX, NY, DEN, ISOCHK

      CALL OLIVER_GET_GRAD_FAST(NOEL, G1X, G1Y, G2X, G2Y,
     &                          G3X, G3Y, HAVE)

      IF (.NOT. HAVE) THEN
         H = CELENT
         IFLAG = 0
         RETURN
      END IF

      ISOCHK = DABS(EXX-EYY) + DABS(GXY)
      IF (ISOCHK .LT. 1.0D-18) THEN
         NX = 1.0D0
         NY = 0.0D0
      ELSE
         THETA = 0.5D0 * DATAN2(GXY, EXX-EYY)
         NX = DCOS(THETA)
         NY = DSIN(THETA)
      END IF

      DEN = DABS(G1X*NX + G1Y*NY) +
     &      DABS(G2X*NX + G2Y*NY) +
     &      DABS(G3X*NX + G3Y*NY)
      IF (DEN .LT. 1.0D-14) DEN = 1.0D-14
      H = 2.0D0 / DEN
      IF (H .LT. 1.0D-12) H = 1.0D-12
      IFLAG = 1

      RETURN
      END

C=======================================================================
C  Fast direct-indexed gradient table reader.
C=======================================================================
      SUBROUTINE OLIVER_GET_GRAD_FAST(NOEL, G1X, G1Y, G2X, G2Y,
     &                                G3X, G3Y, HAVE)

      INCLUDE 'ABA_PARAM.INC'

      INTEGER NOEL
      REAL*8 G1X, G1Y, G2X, G2Y, G3X, G3Y
      LOGICAL HAVE

      INTEGER MAXE, I, IE, IOS, IOS2, NREAD, LU
      PARAMETER (MAXE = 1000000)
      CHARACTER*256 LINE
      LOGICAL INIT, READY, HAS(MAXE), WARNED
      REAL*8 AG1X(MAXE), AG1Y(MAXE), AG2X(MAXE), AG2Y(MAXE)
      REAL*8 AG3X(MAXE), AG3Y(MAXE)
      REAL*8 TG1X, TG1Y, TG2X, TG2Y, TG3X, TG3Y

      SAVE INIT, READY, WARNED, HAS
      SAVE AG1X, AG1Y, AG2X, AG2Y, AG3X, AG3Y
      DATA INIT /.FALSE./
      DATA READY /.FALSE./
      DATA WARNED /.FALSE./

      IF (.NOT. INIT) THEN
         INIT = .TRUE.
         READY = .FALSE.
         WARNED = .FALSE.
         DO I = 1, MAXE
            HAS(I) = .FALSE.
         END DO

         LU = 97
         OPEN(UNIT=LU, FILE='oliver_t3_gradN.dat', STATUS='OLD',
     &        IOSTAT=IOS)
         IF (IOS .NE. 0) THEN
            OPEN(UNIT=LU, FILE='oliver_gradN.dat', STATUS='OLD',
     &           IOSTAT=IOS)
         END IF

         NREAD = 0
         IF (IOS .EQ. 0) THEN
  100       CONTINUE
            READ(LU,'(A)',IOSTAT=IOS) LINE
            IF (IOS .NE. 0) GOTO 200
            IF (LINE .EQ. ' ') GOTO 100
            IF (LINE(1:1) .EQ. '#') GOTO 100
            READ(LINE,*,IOSTAT=IOS2) IE, TG1X, TG1Y, TG2X,
     &           TG2Y, TG3X, TG3Y
            IF (IOS2 .EQ. 0) THEN
               IF (IE .GE. 1 .AND. IE .LE. MAXE) THEN
                  AG1X(IE) = TG1X
                  AG1Y(IE) = TG1Y
                  AG2X(IE) = TG2X
                  AG2Y(IE) = TG2Y
                  AG3X(IE) = TG3X
                  AG3Y(IE) = TG3Y
                  HAS(IE) = .TRUE.
                  NREAD = NREAD + 1
               END IF
            END IF
            GOTO 100
  200       CONTINUE
            CLOSE(LU)
            IF (NREAD .GT. 0) THEN
               READY = .TRUE.
               WRITE(6,*) 'UMAT FAST: Oliver T3 gradN table loaded, n=',
     &              NREAD
            END IF
         END IF

         IF (.NOT. READY) THEN
            WRITE(6,*) 'UMAT WARNING: oliver_t3_gradN.dat not found.'
            WRITE(6,*) 'UMAT WARNING: using Abaqus CELENT fallback.'
         END IF
      END IF

      HAVE = .FALSE.
      IF (READY) THEN
         IF (NOEL .GE. 1 .AND. NOEL .LE. MAXE) THEN
            IF (HAS(NOEL)) THEN
               G1X = AG1X(NOEL)
               G1Y = AG1Y(NOEL)
               G2X = AG2X(NOEL)
               G2Y = AG2Y(NOEL)
               G3X = AG3X(NOEL)
               G3Y = AG3Y(NOEL)
               HAVE = .TRUE.
               RETURN
            END IF
         END IF
      END IF

      G1X = 0.0D0
      G1Y = 0.0D0
      G2X = 0.0D0
      G2Y = 0.0D0
      G3X = 0.0D0
      G3Y = 0.0D0

      IF ((.NOT. HAVE) .AND. READY .AND. (.NOT. WARNED)) THEN
         WRITE(6,*) 'UMAT WARNING: missing element in Oliver table; ',
     &              'CELENT fallback used for missing labels.'
         WARNED = .TRUE.
      END IF

      RETURN
      END
