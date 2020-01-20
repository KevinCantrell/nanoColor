!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been marked up so that f2py can properly wrap    !
! the functions to be used in Python testing.                    !
! Additionally, the "main" function has been removed so that the !
! code can be properly compiled into a shared object.            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c*********** Light Scattering by Spherical Particles *************
c * n-miev3 (n-miev2 + Rayleigh approximation                    *
c *                  + EMT -- n-layered spheres) *               *
c                                                                *
c   Calculations of extinction, scattering, absorption, etc.     *
c  efficiency factors for n-layered spheres.                     *
c................................................................*
c  Input data:          filename: n-miev3.dat                    *
c   n_l              : number of layers (=> 3 !)                 *
c   ri_n(1)   = n-k*i: complex index of refr.(layer 1, innermost)*
c   ri_n(2)   = n-k*i: complex index of refr. (layer 2)          *
c   ..................................................           *
c   ri_n(n_l) = n-k*i: complex index of refr. (layer n_l)        *
c            iv  : relative radius (=0) or volume (=1)           *
c         aa(1)  : relative radius of layer 1 (0 < a1 < 1)       *
c   or    vv(1)  : relative volume of layer 1 (0 < v1 < 1)       *
c         aa(2)  : relative radius of layer 2 (0 < a2 < 1)       *
c   or    vv(2)  : relative volume of layer 2 (0 < v2 < 1)       *
c         ...............................................        *
c         ...............................................        *
c         aa(n_l): relative radius of layer n_l (0 < a(n_l) < 1) *
c   or  vv(n_l-1): relative volume of layer n_l-1                *
c            a   : total radius of a particle (1 + 2 + ...)      *
c!!                vv(n_l) = 1 - [vv(1) + vv(2) + ... vv(n_l-1)] *
c        x_ini   : initial size parameter (2*pi*a/lambda)        *
c           dx   : step                                          *
c        x_fin   : final size parameter                          *
c................................................................*
c  Output data:         filename: n-miev3.out                    *
c   m_eff : effective refractive index                           *
c    Qext : extinction factors                                   *
c    Qsca : scattering factors                                   *
c    Qabs : absorption factors                                   *
c    Qbk  : backscattering factors                               *
c    Qpr  : radiation pressure factors                           *
c   albedo: particle's albedo                                    *
c    g    : asymmetry factor                                     *
c d(a_i)/a: relative radius of layer i (0 < d(a_i)/a < 1)        *
c           (a -- total radius of a particle)                    *
c................................................................*
c NB! In order to treat very large particles,                    *
c     one needs to enlarge the parameter NTERMS.                 *
c................................................................*
c Recursive algorithms of Wu & Wang (Radio Sci. 26, 1393, 1991)  *
c created by N.V. Voshchinnikov                                  *
c (c) 1999 Sobolev Astronomical Institute, St. Petersburg Univ.  *
c*****************************************************************
c

c--------------------------------------------------------------------
c **********   shexqn1 - Spheres: n-layers
c                        Theory: exact
c                        Results: efficiency factors
c March 1999, AI SPbU
c--------------------------------------------------------------------
      SUBROUTINE shexqn1(n_l,ri_n,aa,x,
     *                  qext,qsca,qabs,qbk,qpr,alb,g)
      parameter(n_layers=100, nterms=1000)
      IMPLICIT REAL*8(A-H,O-Q,T-Z), COMPLEX*16(R-S)
cf2py integer,    intent(in)                 :: n_l
cf2py complex*16, intent(in), dimension(n_l) :: ri_n
cf2py real*8,     intent(in), dimension(n_l) :: aa
cf2py real*8,     intent(in)                 :: x
cf2py real*8,     intent(out)                :: qext
cf2py real*8,     intent(out)                :: qsca
cf2py real*8,     intent(out)                :: qabs
cf2py real*8,     intent(out)                :: qbk
cf2py real*8,     intent(out)                :: qpr
cf2py real*8,     intent(out)                :: alb
cf2py real*8,     intent(out)                :: g
      DIMENSION RA(nterms), RB(nterms),
     *          ri_n(n_l), aa(n_l), xx(n_layers),
     *          d1x(nterms), rcx(nterms),
     *          rbb(nterms), rd11(nterms),
     *          rd1(nterms), rd2(nterms), rd3x(nterms),
     *          rrbb(n_layers,nterms),
     *          rrd1(n_layers,nterms), rrd2(n_layers,nterms),
     *          srbb(n_layers,nterms),
     *          srd1(n_layers,nterms), srd2(n_layers,nterms)
      common/fun/ rrbb, rrd1, rrd2, srbb, srd1, srd2,
     *            rd11, rd3x, rcx, d1x

      AX=1.0D0/X
         xx(1) = x * aa(1)
         xx(n_l) = x
        do i = 2, n_l - 1
         xx(i) = 0d0
          do j = 1, i
            xx(i) = xx(i) + aa(j)
          end do
         xx(i) = x * xx(i)
        end do

c d1(x), rd3(x), rc(x)
        NUM = NM(X)
        if(num.gt.nterms) then
         write(*,*) 'nterms, num=', nterms, num
         pause
         stop
        end if
        CALL aax(AX,NUM,d1x)
        CALL cd3x(X,NUM,d1x,rd3x,rcx)

        ari = cdabs(RI_n(1))
        do i = 2, n_l
         ari1 = cdabs(RI_n(i))
         if(ari1.gt.ari) ari = ari1
        end do
        NUM2=NM(ari*X)
        if(num2.gt.nterms) then
         write(*,*) 'nterms, num2=', nterms, num2
         pause
         stop
        end if

c rd11(m_1*x_1)
        if(dimag(ri_n(1))*xx(1).gt.20d0) then
         write(*,*) 'k*x > 20', dimag(ri_n(1))*xx(1)
         pause
        end if
        CALL aa1(ri_n(1)*xx(1),NUM2,rd11)

         do i = 2, n_l
c rd1(m_i*x_i-1), rd2(m_i*x_i-1), rbb(m_i*x_i-1), rcc(m_i*x_i-1),
        if(dimag(ri_n(i))*xx(i-1).gt.20d0) then
         write(*,*) 'k*x > 20', dimag(ri_n(i))*xx(i-1)
         pause
        end if
          CALL bcd(ri_n(i)*xx(i-1),NUM2,rd1,rd2,rbb)
            do j = 1, num2
              rrbb(i,j) = rbb(j)
              rrd1(i,j) = rd1(j)
              rrd2(i,j) = rd2(j)
            end do
c rd1(m_i*x_i), rd2(m_i*x_i), rbb(m_i*x_i), rcc(m_i*x_i),
        if(dimag(ri_n(i))*xx(i).gt.20d0) then
         write(*,*) 'k*x > 20', dimag(ri_n(i))*xx(i)
         pause
        end if
          CALL bcd(ri_n(i)*xx(i),NUM2,rd1,rd2,rbb)
            do j = 1, num2
              srbb(i,j) = rbb(j)
              srd1(i,j) = rd1(j)
              srd2(i,j) = rd2(j)
            end do
         end do

      CALL ABn1(n_l,RI_n,NUM,NUM1,RA,RB)
      CALL QQ1(AX,NUM1,QEXT,QSCA,qbk,qpr,RA,RB)

      qabs=qext-qsca
      alb=qsca/qext
      g=(qext-qpr)/qsca

      RETURN
      END
c--------------------------------------------------------------------
c NM-auxiliary function for AA1 & BESSEL
c    (number NM is calculated using X)
c see: Trudy Astronom. Observ. LGU V.28,P.14,1971
c    for X>1 value of NM was raised
c August 1989, AO LGU
c--------------------------------------------------------------------
      FUNCTION NM(X)
      real*8 x
cf2py real*8,  intent(in)  :: x
cf2py integer, intent(out) :: NM
      IF(X.LT.1) GO TO 11
      IF(X.GT.100) GO TO 12
      NM=1.25*X+15.5
      RETURN
   11 NM=7.5*X+9.0
      RETURN
   12 NM=1.0625*X+28.5
      RETURN
      END
c--------------------------------------------------------------------
c AA1-subroutine for calculations of the ratio of the derivative
c    to the function for Bessel functions of half order with
c    the complex argument: J'(N)/J(N).
c    The calculations are given by the recursive expression
c    ``from top to bottom'' beginning from N=NUM.
c    RU-array of results.
c    A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
c    RI - complex refractive index.
c August 1989, AO LGU
c--------------------------------------------------------------------
      SUBROUTINE AA1(Rx,NUM,RU)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
cf2py complex*16, intent(in)                  :: Rx
cf2py integer,    intent(in)                  :: NUM
cf2py complex*16, dimension(NUM), intent(out) :: RU
      DIMENSION RU(NUM)
      S = 1d0 / rx
      RU(NUM)=(NUM+1.0D0)*S
      NUM1=NUM-1
      DO 13 J=1,NUM1
      I=NUM-J
      I1=I+1
      S1=I1*S
   13 RU(I)=S1-1.0D0/(RU(I1)+S1)
      RETURN
      END
c--------------------------------------------------------------------
c AAx-subroutine for calculations of the ratio of the derivative
c    to the function for Bessel functions of half order with
c    the real argument: J'(N)/J(N).
c    The calculations are given by the recursive expression
c    ``from top to bottom'' beginning from N=NUM.
c    RU-array of results.
c    A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
c March 1999, AI SPbU
c--------------------------------------------------------------------
      SUBROUTINE AAx(A,NUM,RU)
      IMPLICIT REAL*8 (A-H,O-Z)
cf2py real*8,  intent(in)                  :: A
cf2py integer, intent(in)                  :: NUM
cf2py real*8,  dimension(NUM), intent(out) :: RU
      DIMENSION RU(NUM)
      RU(NUM)=(NUM+1.0D0)*a
      NUM1=NUM-1
      DO 13 J=1,NUM1
      I=NUM-J
      I1=I+1
      S1=I1*a
   13 RU(I)=S1-1.0D0/(RU(I1)+S1)
      RETURN
      END
c--------------------------------------------------------------------
c CD3X-subroutine for calculations of the ratio of the derivative
c    to the function for Riccati-Bessel functions of half order with
c    the real argument: zeta'(N)/zeta(N)
c    and the ratio of functions: psi(N)/zeta(N).
c    The calculations are given by the recursive expression
c    ``from bottom to top'' beginning from N=0.
c    rd3x, rcx-arrays of results.
c    X - size parameter
c March 1999, AI SPbU
c--------------------------------------------------------------------
      SUBROUTINE cd3x(x,NUM,d1x,rd3x,rcx)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
cf2py real*8,     intent(in)                  :: x
cf2py integer,    intent(in)                  :: NUM
cf2py real*8,     intent(in),  dimension(NUM) :: d1x
cf2py complex*16, intent(out), dimension(NUM) :: rd3x
cf2py complex*16, intent(out), dimension(NUM) :: rcx
      DIMENSION d1x(NUM), rd3x(NUM), rcx(num)
      S1 = (0d0,1d0)
      ax = 1d0 / x

      rd30 = s1
      rxy = dcos(2d0*x) + s1 * dsin(2d0*x)
      rc0 = -(1d0 - rxy) / (2d0 * rxy)
      rd3x(1) = -ax + 1d0 / (ax - rd30)
      rcx(1) = rc0 * (ax + rd3x(1)) / (ax + d1x(1))

      DO i = 2, NUM
       a1 = I * ax
       rd3x(i) = -a1 + 1d0 / (a1 - rd3x(i-1))
       rcx(i) = rcx(i-1) * (a1 + rd3x(i)) / (a1 + d1x(i))
      end do

      RETURN
      END
c--------------------------------------------------------------------
c BCD-subroutine for calculations of the ratios of the derivative
c    to the function for Riccati-Bessel functions of half order with
c    the complex argument: psi'(N)/psi(N) and khi'(N)/khi(N)
c    and the ratios of functions: psi(N)/khi(N).
c    The calculations are given by the recursive expression
c    ``from bottom to top'' beginning from N=0.
c    rd1, rd2, rbb, rcc-arrays of results.
c    rx - (refr. index) * (size parameter)
c March 1999, AI SPbU
c--------------------------------------------------------------------
      SUBROUTINE bcd(rx,NUM,rd1,rd2,rbb)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      parameter(nterms=1000)
cf2py complex*16, intent(in)                  :: rx
cf2py integer,    intent(in)                  :: NUM
cf2py complex*16, intent(out), dimension(NUM) :: rd1
cf2py complex*16, intent(out), dimension(NUM) :: rd2
cf2py complex*16, intent(out), dimension(NUM) :: rbb
      DIMENSION rd1(NUM), rd2(NUM), rbb(num), rd3(nterms), rcc(nterms)
      S1 = (0d0,1d0)
      x = real(rx)
      y = dimag(rx)
      rx1 = 1d0 / rx

      CALL aa1(rx,NUM,rd1)

c n = 0
      rd30 = s1
      rxy = (dcos(2d0*x) + s1 * dsin(2d0*x))*dexp(-2d0*y)
      rc0 = -(1d0 - rxy) / (2d0 * rxy)
      rb0 = s1 * (1d0 - rxy) / (1d0 + rxy)
c n = 1
      rd3(1) = -rx1 + 1d0 / (rx1 - rd30)
      rcc(1) = rc0 * (rx1 + rd3(1)) / (rx1 + rd1(1))
      rd2(1) = (rcc(1) * rd1(1) - rd3(1)) / (rcc(1) - 1d0)
      rbb(1) = rb0 * (rx1 + rd2(1)) / (rx1 + rd1(1))

      DO i = 2, NUM
       r1 = I * rx1
       rd3(i) = -r1 + 1d0 / (r1 - rd3(i-1))
       rcc(i) = rcc(i-1) * (r1 + rd3(i)) / (r1 + rd1(i))
       rd2(i) = (rcc(i) * rd1(i) - rd3(i)) / (rcc(i) - 1d0)
       rbb(i) = rbb(i-1) * (r1 + rd2(i)) / (r1 + rd1(i))
      end do

      RETURN
      END
c--------------------------------------------------------------------
c ABn1-subroutine for calculations of the complex coefficients
c    A(N), B(N) for n-layered spheres.
c    n_l - number of layers
c    RI_n(i) - complex refractive indices for innermost layer (1),
c    layer2, ... (i = 1, n_l)
c    The coefficients are calculated up to the number NUM1.LE.NUM,
c    for which |A(N)**2+B(N)**2|.LE.10**(-40)
c    RA-array of coefficients A(N), RB-array of coefficients B(N)
c March 1999, AI SPbU
c--------------------------------------------------------------------
      SUBROUTINE ABn1(n_l,RI_n,NUM,NUM1,RA,RB)
      parameter(n_layers=100, nterms=1000)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
cf2py integer,    intent(in)                   :: n_l
cf2py complex*16, intent(in),  dimension(n_l)  :: RI_n
cf2py integer,    intent(in)                   :: NUM
cf2py integer,    intent(out)                  :: NUM1
cf2py complex*16, intent(out), dimension(1000) :: RA
cf2py complex*16, intent(out), dimension(1000) :: RB
      DIMENSION RA(nterms), RB(nterms), ri_n(n_l),
     *          d1x(nterms), rcx(nterms),
     *          rd11(nterms), rd3x(nterms),
     *          rrbb(n_layers,nterms),
     *          rrd1(n_layers,nterms), rrd2(n_layers,nterms),
     *          srbb(n_layers,nterms),
     *          srd1(n_layers,nterms), srd2(n_layers,nterms),
     *          sa(n_layers), sha(n_layers),
     *          sb(n_layers), shb(n_layers)
      common/fun/ rrbb, rrd1, rrd2, srbb, srd1, srd2,
     *            rd11, rd3x, rcx, d1x

c==================
      DO I = 1, NUM

       sa(1) = (0d0,0d0)
       sha(1) = rd11(i)
       sb(1) = (0d0,0d0)
       shb(1) = rd11(i)
c--------
       DO j = 2, n_l
        if(cdabs(ri_n(j)*sha(j-1)-ri_n(j-1)*rrd2(j,i)).eq.0d0) then
        sa(j) = rrbb(j,i) * (ri_n(j) * sha(j-1) - ri_n(j-1) * rrd1(j,i))
     *          / (ri_n(j) * sha(j-1) - ri_n(j-1) * rrd2(j,i) + 1d-30)
        else
        sa(j) = rrbb(j,i) * (ri_n(j) * sha(j-1) - ri_n(j-1) * rrd1(j,i))
     *          / (ri_n(j) * sha(j-1) - ri_n(j-1) * rrd2(j,i))
        end if

        if(cdabs(ri_n(j)*shb(j-1)-ri_n(j-1)*rrd2(j,i)).eq.0d0) then
        sb(j) = rrbb(j,i) * (ri_n(j-1) * shb(j-1) - ri_n(j) * rrd1(j,i))
     *          / (ri_n(j-1) * shb(j-1) - ri_n(j) * rrd2(j,i) + 1d-30)
        else
        sb(j) = rrbb(j,i) * (ri_n(j-1) * shb(j-1) - ri_n(j) * rrd1(j,i))
     *          / (ri_n(j-1) * shb(j-1) - ri_n(j) * rrd2(j,i))
        end if

        if(cdabs(srbb(j,i) - sa(j)).eq.0d0) then
         sha(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sa(j) + 1d-30)
     *           - sa(j) * srd2(j,i) / (srbb(j,i) - sa(j) + 1d-30)
         else
         sha(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sa(j))
     *           - sa(j) * srd2(j,i) / (srbb(j,i) - sa(j))
        end if

        if(cdabs(srbb(j,i) - sb(j)).eq.0d0) then
         shb(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sb(j) + 1d-30)
     *           - sb(j) * srd2(j,i) / (srbb(j,i) - sb(j) + 1d-30)
         else
         shb(j) = srbb(j,i) * srd1(j,i) / (srbb(j,i) - sb(j))
     *           - sb(j) * srd2(j,i) / (srbb(j,i) - sb(j))
        end if

       end do
c--------
* calculations of a(n), b(n)

      RA(I) = rcx(i) * (sha(n_l) - ri_n(n_l) * d1x(i)) /
     *        (sha(n_l) - ri_n(n_l) * rd3x(i))
      RB(I) = rcx(i) * (ri_n(n_l) * shb(n_l) -  d1x(i)) /
     *        (ri_n(n_l) * shb(n_l) - rd3x(i))

      if(cdabs(RA(I))+cdabs(RB(I)).LE.1D-40) GO TO 12
      end do
   12 NUM1=I
      RETURN
      END
c--------------------------------------------------------------------
c QQ1-subroutine for calculations of the efficiency factors for
c     extinction (QEXT), scattering (QSCA), backscattering (QBK)
c     and radiation pressure (QPR) for spherical particles.
c August 1989, AO LGU
c--------------------------------------------------------------------
      SUBROUTINE QQ1(A,NUM,QEXT,QSCA,qbk,qpr,RA,RB)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
cf2py real*8,     intent(in)                 :: A
cf2py integer,    intent(in)                 :: NUM
cf2py real*8,     intent(out)                :: QEXT
cf2py real*8,     intent(out)                :: QSCA
cf2py real*8,     intent(out)                :: qbk
cf2py real*8,     intent(out)                :: qpr
cf2py complex*16, intent(in), dimension(NUM) :: RA
cf2py complex*16, intent(in), dimension(NUM) :: RB
      DIMENSION RA(NUM),RB(NUM)
      B=2.0D0*A*A
      C=0.0D0
      D=0.0D0
      s=(0d0,0d0)
      r=(0d0,0d0)
      N=1
      DO 11 I=1,NUM-1
      N=N+2
      r=r+(i+0.5d0)*(-1)**i*(ra(i)-rb(i))
      s=s+i*(i+2d0)/(i+1d0)*(ra(i)*dconjg(ra(i+1))
     *  +rb(i)*dconjg(rb(i+1)))+n/i/(i+1d0)*(ra(i)*dconjg(rb(i)))
      C=C+N*(RA(I)+RB(I))
   11 D=D+N*(RA(I)*DCONJG(RA(I))+RB(I)*DCONJG(RB(I)))
      QEXT=B*C
      QSCA=B*D
      qbk=2d0*b*r*dconjg(r)
      qpr=qext-2d0*b*s
      RETURN
      END
