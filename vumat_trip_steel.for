c=====================================================================
c  VUMAT for TRIP Steel - Abaqus/Explicit
c  ---------------------------------------------------------------
c  Coupled large-deformation model (homogenised material point):
c    1. Isotropic Elasticity       (Hooke 3D, co-rotational/hypoelastic)
c    2. Crystal Plasticity FCC     (12 slip systems {111}<110>, Voce)
c    3. Martensite Transformation  (Olson-Cohen + JMAK, vol. eigenstrain)
c    4. Twinning                   (12 sys {111}<112>, CRSS + TANH)
c    5. Twin <-> Martensite coupling
c  ---------------------------------------------------------------
c  (Hill48 / Barlat macroscopic yield criteria are NOT included in this
c   version: in a crystal-plasticity framework the per-slip-system CRSS
c   plays the role of the yield surface via the Schmid + power law.)
c  ---------------------------------------------------------------
c  MODELLING NOTES
c  *  Yielding/flow is governed PER SLIP SYSTEM by the resolved shear
c     vs the CRSS (Schmid + power law) - there is no macroscopic yield
c     surface in this model.
c  *  This is a HOMOGENISED single-point model. Austenite (FCC) carries
c     the crystal plasticity. Martensite (BCT) is represented by its
c     volume fraction f_M and its transformation (volumetric) eigen-
c     strain + hardening contribution - NOT by a separate set of BCT
c     slip systems. A full two-phase (FCC+BCT) scheme is a possible
c     extension.
c  *  Temperature is held fixed (default 298 K). Olson-Cohen / JMAK
c     temperature dependence can be folded into props(12) externally.
c  *  Twinning is UNIDIRECTIONAL: a twin system activates only when its
c     resolved shear tau^t > 0 AND tau^t > tau_twin. Verify the sign
c     (sense) of each twin direction against your own crystallography.
c  ---------------------------------------------------------------
c  State Variables (nstatev = 30):
c    sv(1  :12) = gamma_s       slip shear strain per system
c    sv(13 :24) = tauCR_s       CRSS per slip system
c    sv(25)     = f_M           martensite volume fraction
c    sv(26)     = f_tw          twin volume fraction
c    sv(27)     = eps_p_eq      equivalent plastic strain (slip+twin)
c    sv(28)     = Nsb           shear band density (Olson-Cohen)
c    sv(29)     = eps_M         accumulated transformation strain
c    sv(30)     = T             temperature (K)
c  ---------------------------------------------------------------
c  Material Properties (nprops = 20):
c    props(1)  = E              Young modulus (MPa)
c    props(2)  = nu             Poisson ratio
c    props(3)  = gam0           reference shear strain rate (1/s)
c    props(4)  = mexp           rate sensitivity exponent m (0<m<1)
c    props(5)  = tau0           initial CRSS (MPa)
c    props(6)  = taus           saturation CRSS (MPa)
c    props(7)  = h0             hardening modulus (MPa)
c    props(8)  = alpha_OC       Olson-Cohen alpha (shear-band rate)
c    props(9)  = beta_OC        Olson-Cohen beta  (nucleation)
c    props(10) = n_OC           Olson-Cohen n exponent
c    props(11) = eps0_OC        Olson-Cohen reference plastic strain
c    props(12) = K_J            JMAK rate constant (T-corrected ext.)
c    props(13) = n_J            JMAK exponent
c    props(14) = T_ref          JMAK reference temperature (K)
c    props(15) = dVoV           volumetric transf. strain dV/V (~0.02-0.04)
c    props(16) = tau_twin       critical twin shear stress tau_c^tw (MPa)
c    props(17) = gam_twin       characteristic twin shear (~0.707 FCC)
c    props(18) = A_tw           TANH saturation amplitude (max f_tw, <=1)
c    props(19) = B_tw           TANH transition width Dtau (MPa)
c    props(20) = alpha_TM       twin->martensite coupling coeff. (0=off)
c=====================================================================

      subroutine vumat(
     1 nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2 stepTime, totalTime, dt, cmname, coordMp, charLength,
     3 props, density, strainInc, relSpinInc,
     4 tempOld, stretchOld, defgradOld, fieldOld,
     5 stressOld, stateOld, enerInternOld, enerInelasOld,
     6 tempNew, stretchNew, defgradNew, fieldNew,
     7 stressNew, stateNew, enerInternNew, enerInelasNew)

      include 'vaba_param.inc'

c----- Abaqus VUMAT arguments -----------------------------------------
      character*80 cmname
      integer nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal
      double precision stepTime, totalTime, dt
      double precision coordMp(3,nblock), charLength(nblock)
      double precision props(nprops), density(nblock)
      double precision strainInc(ndir+nshr,nblock)
      double precision relSpinInc(nshr,nblock)
      double precision tempOld(nblock), tempNew(nblock)
      double precision stretchOld(3,3,nblock), stretchNew(3,3,nblock)
      double precision defgradOld(3,3,nblock), defgradNew(3,3,nblock)
      double precision fieldOld(nfieldv,nblock),fieldNew(nfieldv,nblock)
      double precision stressOld(ndir+nshr,nblock)
      double precision stressNew(ndir+nshr,nblock)
      double precision stateOld(nstatev,nblock)
      double precision stateNew(nstatev,nblock)
      double precision enerInternOld(nblock), enerInternNew(nblock)
      double precision enerInelasOld(nblock), enerInelasNew(nblock)

c----- Constants ------------------------------------------------------
      integer    NSLIP, NTWIN, NSTATEV_MIN
      parameter (NSLIP = 12)
      parameter (NTWIN = 12)
      parameter (NSTATEV_MIN = 30)

      double precision zero, one, two, three, half, small
      parameter (zero=0.d0, one=1.d0, two=2.d0,
     &           three=3.d0, half=0.5d0, small=1.d-12)

c  Numerical safety caps (explicit stability)
      double precision PWMAX, DGMAX, DFTWMAX, DFMMAX
      parameter (PWMAX = 1.d8)
      parameter (DGMAX = 0.01d0)
c  Max volume-fraction increase per increment (transformation/twinning)
c  spreads rapid transformation over several increments -> robust to
c  aggressive A_tw / B_tw and avoids single-increment strain spikes.
      parameter (DFTWMAX = 0.02d0)
      parameter (DFMMAX  = 0.02d0)

c----- Material parameters --------------------------------------------
      double precision E, nu, gam0, mexp, tau0, taus, h0
      double precision alpha_OC, beta_OC, n_OC, eps0_OC
      double precision K_J, n_J, T_ref
      double precision dVoV
      double precision tau_twin, gam_twin, A_tw, B_tw, alpha_TM
      double precision lambda, mu

c----- Elastic stiffness D(6x6) - Isotropic Hooke --------------------
      double precision D(6,6)

c----- Slip system variables ------------------------------------------
      double precision Ms(3,3,NSLIP)
      double precision Tau_s(NSLIP)
      double precision Gamma_old(NSLIP), Gamma_new(NSLIP)
      double precision TauCR_old(NSLIP), TauCR_new(NSLIP)
      double precision Dgamma(NSLIP)

c----- Twin system variables ------------------------------------------
      double precision Mtw(3,3,NTWIN)
      double precision Tau_tw(NTWIN), wtw(NTWIN)

c----- Stress / strain arrays -----------------------------------------
      double precision de(6), ds(6), de_el(6)
      double precision Einc(3,3), Ep_cp(3,3), Ep_tw(3,3), Ep_tot(3,3)
      double precision Eel(3,3)
      double precision SigmaOld(3,3), SigmaTrial(3,3)

c----- Phase transformation variables ---------------------------------
      double precision f_M_old, f_M_new, f_M_base
      double precision f_tw_old, f_tw_new, f_tw_target, f_tw_prelim
      double precision eps_p_eq_old, eps_p_eq_new
      double precision Nsb_old, Nsb_new
      double precision eps_M_old, eps_M_new
      double precision T_loc

c----- Martensite / twin work scalars ---------------------------------
      double precision f_M_OC, f_M_JMAK
      double precision dNsb, df_M, dftw, dftw_eff
      double precision deps_M_vol, Bulk
      double precision tau_tw_max, wsum, twfac

c----- Loop indices ---------------------------------------------------
      integer ib, i, j, islip, itw, ntens, ic

c----- Temporary scalars ----------------------------------------------
      double precision at, r, signTau, gdot, pw
      double precision arg, fac, gamma_eq
      double precision stress6(6)
      double precision deps_p_eq, Ep_norm, cap_tw
      double precision sfac, dev_inc, dev_p, trace_e

c=====================================================================
c  Read material properties
c=====================================================================
      E         = props(1)
      nu        = props(2)
      gam0      = props(3)
      mexp      = props(4)
      tau0      = props(5)
      taus      = props(6)
      h0        = props(7)
      alpha_OC  = props(8)
      beta_OC   = props(9)
      n_OC      = props(10)
      eps0_OC   = props(11)
      K_J       = props(12)
      n_J       = props(13)
      T_ref     = props(14)
      dVoV      = props(15)
      tau_twin  = props(16)
      gam_twin  = props(17)
      A_tw      = props(18)
      B_tw      = props(19)
      alpha_TM  = props(20)

c----- Lame constants -------------------------------------------------
      mu     = E / (two*(one + nu))
      lambda = E*nu / ((one + nu)*(one - two*nu))

c  Bulk modulus: K = lambda + 2/3*mu
      Bulk   = lambda + two*mu/three

c=====================================================================
c  Build Isotropic Elastic Stiffness D(6x6) - Hooke 3D
c  Voigt order: [s11, s22, s33, s12, s13, s23]
c  NOTE: Abaqus VUMAT strainInc shear terms are TENSORIAL (eps_12),
c        so the shear stiffness must be 2*mu  (sigma_12 = 2*mu*eps_12).
c=====================================================================
      do i = 1,6
        do j = 1,6
          D(i,j) = zero
        end do
      end do
      D(1,1) = lambda + two*mu
      D(2,2) = lambda + two*mu
      D(3,3) = lambda + two*mu
      D(1,2) = lambda
      D(2,1) = lambda
      D(1,3) = lambda
      D(3,1) = lambda
      D(2,3) = lambda
      D(3,2) = lambda
      D(4,4) = two*mu
      D(5,5) = two*mu
      D(6,6) = two*mu

c----- FCC Schmid tensors (slip + twin) -------------------------------
      call init_fcc_schmid(Ms)
      call init_fcc_twin(Mtw)

      ntens = ndir + nshr

c=====================================================================
c  Main block loop
c=====================================================================
      do ib = 1, nblock

        tempNew(ib) = tempOld(ib)
        T_loc = tempOld(ib)
        if (T_loc .le. small) T_loc = 298.d0

        do i = 1, nfieldv
          fieldNew(i,ib) = fieldOld(i,ib)
        end do

        do i = 1,6
          if (i .le. ntens) then
            de(i) = strainInc(i,ib)
          else
            de(i) = zero
          end if
        end do

c=====================================================================
c  CASE 1: stepTime = 0 -> purely elastic initialisation
c=====================================================================
        if (stepTime .eq. zero) then

          do i = 1,6
            ds(i) = zero
            do j = 1,6
              ds(i) = ds(i) + D(i,j)*de(j)
            end do
          end do

          do i = 1, ntens
            stressNew(i,ib) = stressOld(i,ib) + ds(i)
          end do

          do islip = 1, NSLIP
            Gamma_old(islip) = stateOld(islip,ib)
            TauCR_old(islip) = stateOld(NSLIP+islip,ib)
            if (TauCR_old(islip) .le. zero) TauCR_old(islip) = tau0
            stateNew(islip,ib)       = Gamma_old(islip)
            stateNew(NSLIP+islip,ib) = TauCR_old(islip)
          end do
          stateNew(25,ib) = stateOld(25,ib)
          stateNew(26,ib) = stateOld(26,ib)
          stateNew(27,ib) = stateOld(27,ib)
          stateNew(28,ib) = stateOld(28,ib)
          stateNew(29,ib) = stateOld(29,ib)
          stateNew(30,ib) = T_loc

          enerInternNew(ib) = enerInternOld(ib)
          enerInelasNew(ib) = enerInelasOld(ib)

          cycle
        end if

c=====================================================================
c  CASE 2: stepTime > 0 -> full coupled model
c=====================================================================

        if (totalTime .eq. zero) then
          do islip = 1, NSLIP
            stateNew(islip,ib)       = zero
            stateNew(NSLIP+islip,ib) = tau0
          end do
          stateNew(25,ib) = zero
          stateNew(26,ib) = zero
          stateNew(27,ib) = zero
          stateNew(28,ib) = zero
          stateNew(29,ib) = zero
          stateNew(30,ib) = T_loc
        end if

c  Read old state
        do islip = 1, NSLIP
          Gamma_old(islip) = stateOld(islip,ib)
          TauCR_old(islip) = stateOld(NSLIP+islip,ib)
          if (TauCR_old(islip) .le. zero) TauCR_old(islip) = tau0
        end do
        f_M_old      = stateOld(25,ib)
        f_tw_old     = stateOld(26,ib)
        eps_p_eq_old = stateOld(27,ib)
        Nsb_old      = stateOld(28,ib)
        eps_M_old    = stateOld(29,ib)

        if (f_M_old  .lt. zero) f_M_old  = zero
        if (f_M_old  .gt. one)  f_M_old  = one
        if (f_tw_old .lt. zero) f_tw_old = zero
        if (f_tw_old .gt. one)  f_tw_old = one

c---------------------------------------------------------------------
c  STEP A: Elastic predictor (Hooke 3D)
c  sigma_trial = sigma_old + D : de
c---------------------------------------------------------------------
        do i = 1,6
          ds(i) = zero
          do j = 1,6
            ds(i) = ds(i) + D(i,j)*de(j)
          end do
        end do

        do i = 1, ntens
          stress6(i) = stressOld(i,ib) + ds(i)
        end do
        do i = ntens+1, 6
          stress6(i) = zero
        end do

        call mat_from_vec(stress6, SigmaTrial)
        call mat_from_vec(stressOld(1,ib), SigmaOld)

c---------------------------------------------------------------------
c  STEP B: Crystal Plasticity FCC (12 slip systems)
c  ---------------------------------------------------------------
c    Power-law:  gdot_s = gam0*|tau_s/tauCR_s|^(1/m)*sign(tau_s)
c    Voce:       tauCR  = tau0+(taus-tau0)*(1-exp(-h0*Gam/(taus-tau0)))
c    Ep_cp      = sum_s Ms_s * Dgamma_s
c---------------------------------------------------------------------
        do islip = 1, NSLIP
          Tau_s(islip) = zero
          do i = 1,3
            do j = 1,3
              Tau_s(islip) = Tau_s(islip)
     &                     + Ms(i,j,islip)*SigmaTrial(i,j)
            end do
          end do
        end do

        do islip = 1, NSLIP
          at = dabs(Tau_s(islip))
          if (at .le. small) then
            Dgamma(islip) = zero
          else
            if (TauCR_old(islip) .le. small)
     &          TauCR_old(islip) = tau0
            r = at / TauCR_old(islip)
            if (r .lt. small) r = small
            signTau = Tau_s(islip) / at
            pw = r**(one/mexp)
            if (pw .gt. PWMAX) pw = PWMAX
            gdot = gam0 * pw * signTau
            Dgamma(islip) = gdot * dt
c           Cap increment to keep explicit integration stable
            if (Dgamma(islip) .gt.  DGMAX) Dgamma(islip) =  DGMAX
            if (Dgamma(islip) .lt. -DGMAX) Dgamma(islip) = -DGMAX
          end if
        end do

        call zero_tensor(Ep_cp)
        do islip = 1, NSLIP
          do i = 1,3
            do j = 1,3
              Ep_cp(i,j) = Ep_cp(i,j)
     &                   + Ms(i,j,islip)*Dgamma(islip)
            end do
          end do
        end do

c---------------------------------------------------------------------
c  GLOBAL ANTI-OVERSHOOT LIMITER (explicit stability)
c  ---------------------------------------------------------------
c  Forward-Euler slip can predict more plastic strain than was imposed
c  in the increment -> stress reverses sign -> element collapse.
c  Cap the (deviatoric) plastic strain increment so it never exceeds
c  the imposed deviatoric strain increment magnitude. This is the
c  rate-independent fully-plastic limit and removes the instability.
c---------------------------------------------------------------------
        call mat_from_vec(de, Einc)
        trace_e = (Einc(1,1)+Einc(2,2)+Einc(3,3))/three
        dev_inc = zero
        do i = 1,3
          do j = 1,3
            if (i .eq. j) then
              dev_inc = dev_inc + (Einc(i,j)-trace_e)**2
            else
              dev_inc = dev_inc + Einc(i,j)**2
            end if
          end do
        end do
        dev_inc = dsqrt(dev_inc)

        dev_p = zero
        do i = 1,3
          do j = 1,3
            dev_p = dev_p + Ep_cp(i,j)**2
          end do
        end do
        dev_p = dsqrt(dev_p)

        if (dev_p .gt. dev_inc .and. dev_p .gt. small) then
          sfac = dev_inc / dev_p
          do islip = 1, NSLIP
            Dgamma(islip) = Dgamma(islip)*sfac
          end do
          do i = 1,3
            do j = 1,3
              Ep_cp(i,j) = Ep_cp(i,j)*sfac
            end do
          end do
        end if

c  Voce hardening update (per system)
        do islip = 1, NSLIP
          Gamma_new(islip) = Gamma_old(islip) + dabs(Dgamma(islip))
          gamma_eq = Gamma_new(islip)
          if (taus .gt. tau0 + small) then
            arg = -h0 * gamma_eq / (taus - tau0)
            if (arg .lt. -50.d0) arg = -50.d0
            fac = one - dexp(arg)
            TauCR_new(islip) = tau0 + (taus - tau0)*fac
          else
            TauCR_new(islip) = tau0 + h0*gamma_eq
          end if
        end do

c---------------------------------------------------------------------
c  STEP D: Twinning (CRSS criterion + TANH kinetics)  [pseudo-slip]
c  ---------------------------------------------------------------
c    Resolved twin shear:  tau^t = Mtw_t : sigma_trial
c    Unidirectional: only tau^t > 0 (and > tau_twin) activates.
c    TANH volume-fraction target (driven by the most stressed system):
c      f_tw_target = A_tw * 0.5*(1 + tanh((tau_max - tau_twin)/B_tw))
c    Irreversible growth, capped by available austenite (1 - f_M).
c    Twin plastic strain (distributed over active systems by weight):
c      Ep_tw = sum_t (w_t/W) * dftw * gam_twin * Mtw_t
c---------------------------------------------------------------------
        do itw = 1, NTWIN
          Tau_tw(itw) = zero
          do i = 1,3
            do j = 1,3
              Tau_tw(itw) = Tau_tw(itw)
     &                    + Mtw(i,j,itw)*SigmaTrial(i,j)
            end do
          end do
        end do

c  Active-system weights and maximum positive resolved shear
        tau_tw_max = zero
        wsum       = zero
        do itw = 1, NTWIN
          if (Tau_tw(itw) .gt. tau_twin) then
            wtw(itw) = Tau_tw(itw) - tau_twin
          else
            wtw(itw) = zero
          end if
          wsum = wsum + wtw(itw)
          if (Tau_tw(itw) .gt. tau_tw_max) tau_tw_max = Tau_tw(itw)
        end do

c  TANH target twin fraction
        if (B_tw .gt. small) then
          arg = (tau_tw_max - tau_twin) / B_tw
          if (arg .gt.  50.d0) arg =  50.d0
          if (arg .lt. -50.d0) arg = -50.d0
          f_tw_target = A_tw * half*(one + dtanh(arg))
        else
          if (tau_tw_max .gt. tau_twin) then
            f_tw_target = A_tw
          else
            f_tw_target = zero
          end if
        end if

c  Cap by available austenite (use old f_M to avoid circularity)
        cap_tw = one - f_M_old
        if (cap_tw .lt. zero) cap_tw = zero
        if (f_tw_target .gt. cap_tw) f_tw_target = cap_tw

c  Irreversible growth
        dftw = f_tw_target - f_tw_old
        if (dftw .lt. zero) dftw = zero
        if (dftw .gt. DFTWMAX) dftw = DFTWMAX
        f_tw_prelim = f_tw_old + dftw

c---------------------------------------------------------------------
c  STEP C: Martensite Transformation (Olson-Cohen + JMAK)
c  ---------------------------------------------------------------
c  C1. Olson-Cohen shear-band nucleation:
c        Nsb_new = Nsb_old + alpha_OC*(1-Nsb_old)*deps_p_eq/eps0_OC
c        f_M_OC  = 1 - exp(-beta_OC * Nsb_new^n_OC)
c  C2. JMAK kinetics:
c        f_M_JMAK = 1 - exp(-K_J * eps_p_eq_new^n_J)
c  C3. Combine as independent mechanisms, then add twin coupling:
c        f_M_base = 1 - (1-f_M_OC)*(1-f_M_JMAK)
c        f_M_new  = f_M_base + alpha_TM*dftw*(1 - f_M_base)
c  C4. Exact volumetric eigenstrain: dsigma_ii = -K*df_M*dVoV
c---------------------------------------------------------------------
c  Equivalent plastic strain increment (slip only drives nucleation)
        Ep_norm = zero
        do i = 1,3
          do j = 1,3
            Ep_norm = Ep_norm + Ep_cp(i,j)*Ep_cp(i,j)
          end do
        end do
        deps_p_eq    = dsqrt(two/three * Ep_norm)
        eps_p_eq_new = eps_p_eq_old + deps_p_eq

c  --- C1: Olson-Cohen shear band evolution ---
        if (eps0_OC .gt. small) then
          dNsb = alpha_OC * (one - Nsb_old)
     &         * (deps_p_eq / eps0_OC)
        else
          dNsb = zero
        end if
        Nsb_new = Nsb_old + dNsb
        if (Nsb_new .gt. one)  Nsb_new = one
        if (Nsb_new .lt. zero) Nsb_new = zero

        if (Nsb_new .gt. small .and. beta_OC .gt. small) then
          arg = -beta_OC * (Nsb_new**n_OC)
          if (arg .lt. -50.d0) arg = -50.d0
          f_M_OC = one - dexp(arg)
        else
          f_M_OC = zero
        end if

c  --- C2: JMAK kinetics ---
        if (eps_p_eq_new .gt. small .and. K_J .gt. small) then
          arg = -K_J * (eps_p_eq_new**n_J)
          if (arg .lt. -50.d0) arg = -50.d0
          f_M_JMAK = one - dexp(arg)
        else
          f_M_JMAK = zero
        end if

c  --- C3: Combined base fraction (independent events) ---
        f_M_base = one - (one - f_M_OC)*(one - f_M_JMAK)

c  Twin -> martensite coupling: extra nucleation at twin boundaries
        f_M_base = f_M_base + alpha_TM*dftw*(one - f_M_base)

c  Monotonicity: martensite cannot revert
        if (f_M_base .gt. f_M_old) then
          f_M_new = f_M_base
        else
          f_M_new = f_M_old
        end if
        if (f_M_new .gt. one) f_M_new = one
c  Rate cap: limit martensite growth per increment (explicit stability)
        if (f_M_new - f_M_old .gt. DFMMAX) f_M_new = f_M_old + DFMMAX
        df_M = f_M_new - f_M_old

c  Re-cap twin so that f_tw + f_M <= 1, rescale twin strain accordingly
        cap_tw = one - f_M_new
        if (cap_tw .lt. zero) cap_tw = zero
        f_tw_new = f_tw_prelim
        if (f_tw_new .gt. cap_tw) f_tw_new = cap_tw
        dftw_eff = f_tw_new - f_tw_old
        if (dftw_eff .lt. zero) dftw_eff = zero

c  --- Twin plastic strain tensor (distributed over active systems) ---
        call zero_tensor(Ep_tw)
        if (wsum .gt. small .and. dftw_eff .gt. small) then
          do itw = 1, NTWIN
            if (wtw(itw) .gt. zero) then
              twfac = (wtw(itw)/wsum) * dftw_eff * gam_twin
              do i = 1,3
                do j = 1,3
                  Ep_tw(i,j) = Ep_tw(i,j) + twfac*Mtw(i,j,itw)
                end do
              end do
            end if
          end do
        end if

c---------------------------------------------------------------------
c  STEP E: Assemble inelastic strain and update stress
c  ---------------------------------------------------------------
c    Ep_tot = Ep_cp + Ep_tw
c    sigma  = sigma_old + D:(de - Ep_tot)
c    then subtract martensite volumetric eigenstrain
c---------------------------------------------------------------------
        do i = 1,3
          do j = 1,3
            Ep_tot(i,j) = Ep_cp(i,j) + Ep_tw(i,j)
          end do
        end do

        do i = 1,3
          do j = 1,3
            Eel(i,j) = Einc(i,j) - Ep_tot(i,j)
          end do
        end do

        call vec_from_mat(Eel, de_el)
        do i = 1,6
          ds(i) = zero
          do j = 1,6
            ds(i) = ds(i) + D(i,j)*de_el(j)
          end do
        end do

c  Martensite volumetric eigenstrain correction (isotropic)
c    eps_tr_vol(per dir) = df_M*dVoV/3 ; dsigma_ii = -3K*eps_tr_vol
        deps_M_vol = df_M * dVoV / three
        ds(1) = ds(1) - three * Bulk * deps_M_vol
        ds(2) = ds(2) - three * Bulk * deps_M_vol
        ds(3) = ds(3) - three * Bulk * deps_M_vol

        eps_M_new = eps_M_old + dabs(deps_M_vol)*three

        do i = 1, ntens
          stressNew(i,ib) = stressOld(i,ib) + ds(i)
        end do

c---------------------------------------------------------------------
c  STEP F: Write state variables
c---------------------------------------------------------------------
        do islip = 1, NSLIP
          stateNew(islip,ib)       = Gamma_new(islip)
          stateNew(NSLIP+islip,ib) = TauCR_new(islip)
        end do
        stateNew(25,ib) = f_M_new
        stateNew(26,ib) = f_tw_new
        stateNew(27,ib) = eps_p_eq_new
        stateNew(28,ib) = Nsb_new
        stateNew(29,ib) = eps_M_new
        stateNew(30,ib) = T_loc

c  Internal / inelastic energy (approx.: slip + twin dissipation)
        enerInternNew(ib) = enerInternOld(ib)
        enerInelasNew(ib) = enerInelasOld(ib) + deps_p_eq

      end do
c  End block loop

      return
      end

c=====================================================================
c  HELPER: FCC slip systems {111}<110> - 12 systems
c=====================================================================
      subroutine init_fcc_schmid(Ms)
      include 'vaba_param.inc'
      integer NSLIP
      parameter (NSLIP = 12)
      double precision Ms(3,3,NSLIP)
      double precision n(3), s(3)
      double precision sqrt2, sqrt3
      integer i,j,islip

      sqrt2 = dsqrt(2.d0)
      sqrt3 = dsqrt(3.d0)

      do islip = 1, NSLIP
        do i = 1,3
          do j = 1,3
            Ms(i,j,islip) = 0.d0
          end do
        end do
      end do

c  Plane (1 1 1)
      n(1)= 1.d0/sqrt3
      n(2)= 1.d0/sqrt3
      n(3)= 1.d0/sqrt3
      s(1)= 0.d0
      s(2)=-1.d0/sqrt2
      s(3)= 1.d0/sqrt2
      call build_schmid(n,s,Ms,1)
      s(1)= 1.d0/sqrt2
      s(2)= 0.d0
      s(3)=-1.d0/sqrt2
      call build_schmid(n,s,Ms,2)
      s(1)=-1.d0/sqrt2
      s(2)= 1.d0/sqrt2
      s(3)= 0.d0
      call build_schmid(n,s,Ms,3)

c  Plane (1 -1 1)
      n(1)= 1.d0/sqrt3
      n(2)=-1.d0/sqrt3
      n(3)= 1.d0/sqrt3
      s(1)= 0.d0
      s(2)= 1.d0/sqrt2
      s(3)= 1.d0/sqrt2
      call build_schmid(n,s,Ms,4)
      s(1)= 1.d0/sqrt2
      s(2)= 0.d0
      s(3)= 1.d0/sqrt2
      call build_schmid(n,s,Ms,5)
      s(1)=-1.d0/sqrt2
      s(2)=-1.d0/sqrt2
      s(3)= 0.d0
      call build_schmid(n,s,Ms,6)

c  Plane (-1 1 1)
      n(1)=-1.d0/sqrt3
      n(2)= 1.d0/sqrt3
      n(3)= 1.d0/sqrt3
      s(1)= 0.d0
      s(2)= 1.d0/sqrt2
      s(3)=-1.d0/sqrt2
      call build_schmid(n,s,Ms,7)
      s(1)= 1.d0/sqrt2
      s(2)= 0.d0
      s(3)= 1.d0/sqrt2
      call build_schmid(n,s,Ms,8)
      s(1)=-1.d0/sqrt2
      s(2)= 1.d0/sqrt2
      s(3)= 0.d0
      call build_schmid(n,s,Ms,9)

c  Plane (1 1 -1)
      n(1)= 1.d0/sqrt3
      n(2)= 1.d0/sqrt3
      n(3)=-1.d0/sqrt3
      s(1)= 0.d0
      s(2)= 1.d0/sqrt2
      s(3)= 1.d0/sqrt2
      call build_schmid(n,s,Ms,10)
      s(1)= 1.d0/sqrt2
      s(2)= 0.d0
      s(3)=-1.d0/sqrt2
      call build_schmid(n,s,Ms,11)
      s(1)=-1.d0/sqrt2
      s(2)=-1.d0/sqrt2
      s(3)= 0.d0
      call build_schmid(n,s,Ms,12)

      return
      end

c=====================================================================
c  HELPER: FCC twin systems {111}<112> - 12 systems
c  ---------------------------------------------------------------
c  Twinning is UNIDIRECTIONAL. The shear direction sign below sets the
c  twinning sense; verify against your crystallography/texture if you
c  need the exact polar sense. n.s = 0 holds for every pair.
c=====================================================================
      subroutine init_fcc_twin(Mtw)
      include 'vaba_param.inc'
      integer NTWIN
      parameter (NTWIN = 12)
      double precision Mtw(3,3,NTWIN)
      double precision n(3), s(3)
      double precision sqrt3, sqrt6
      integer i,j,itw

      sqrt3 = dsqrt(3.d0)
      sqrt6 = dsqrt(6.d0)

      do itw = 1, NTWIN
        do i = 1,3
          do j = 1,3
            Mtw(i,j,itw) = 0.d0
          end do
        end do
      end do

c  Plane (1 1 1)
      n(1)= 1.d0/sqrt3
      n(2)= 1.d0/sqrt3
      n(3)= 1.d0/sqrt3
      s(1)= 1.d0/sqrt6
      s(2)= 1.d0/sqrt6
      s(3)=-2.d0/sqrt6
      call build_schmid(n,s,Mtw,1)
      s(1)= 1.d0/sqrt6
      s(2)=-2.d0/sqrt6
      s(3)= 1.d0/sqrt6
      call build_schmid(n,s,Mtw,2)
      s(1)=-2.d0/sqrt6
      s(2)= 1.d0/sqrt6
      s(3)= 1.d0/sqrt6
      call build_schmid(n,s,Mtw,3)

c  Plane (1 -1 1)
      n(1)= 1.d0/sqrt3
      n(2)=-1.d0/sqrt3
      n(3)= 1.d0/sqrt3
      s(1)= 2.d0/sqrt6
      s(2)= 1.d0/sqrt6
      s(3)=-1.d0/sqrt6
      call build_schmid(n,s,Mtw,4)
      s(1)= 1.d0/sqrt6
      s(2)= 2.d0/sqrt6
      s(3)= 1.d0/sqrt6
      call build_schmid(n,s,Mtw,5)
      s(1)=-1.d0/sqrt6
      s(2)= 1.d0/sqrt6
      s(3)= 2.d0/sqrt6
      call build_schmid(n,s,Mtw,6)

c  Plane (-1 1 1)
      n(1)=-1.d0/sqrt3
      n(2)= 1.d0/sqrt3
      n(3)= 1.d0/sqrt3
      s(1)= 1.d0/sqrt6
      s(2)= 2.d0/sqrt6
      s(3)=-1.d0/sqrt6
      call build_schmid(n,s,Mtw,7)
      s(1)= 2.d0/sqrt6
      s(2)= 1.d0/sqrt6
      s(3)= 1.d0/sqrt6
      call build_schmid(n,s,Mtw,8)
      s(1)= 1.d0/sqrt6
      s(2)=-1.d0/sqrt6
      s(3)= 2.d0/sqrt6
      call build_schmid(n,s,Mtw,9)

c  Plane (1 1 -1)
      n(1)= 1.d0/sqrt3
      n(2)= 1.d0/sqrt3
      n(3)=-1.d0/sqrt3
      s(1)= 1.d0/sqrt6
      s(2)= 1.d0/sqrt6
      s(3)= 2.d0/sqrt6
      call build_schmid(n,s,Mtw,10)
      s(1)= 2.d0/sqrt6
      s(2)=-1.d0/sqrt6
      s(3)= 1.d0/sqrt6
      call build_schmid(n,s,Mtw,11)
      s(1)=-1.d0/sqrt6
      s(2)= 2.d0/sqrt6
      s(3)= 1.d0/sqrt6
      call build_schmid(n,s,Mtw,12)

      return
      end

c=====================================================================
c  HELPER: Build symmetric Schmid tensor Ms = 0.5*(s x n + n x s)
c=====================================================================
      subroutine build_schmid(n,s,Ms,idx)
      include 'vaba_param.inc'
      double precision n(3), s(3), Ms(3,3,*)
      integer idx, i, j
      do i = 1,3
        do j = 1,3
          Ms(i,j,idx) = 0.5d0*(s(i)*n(j) + n(i)*s(j))
        end do
      end do
      return
      end

c=====================================================================
c  HELPER: Voigt vector (6) -> symmetric 3x3 tensor
c  VUMAT order: [11, 22, 33, 12, 23, 13]  (shear terms TENSORIAL)
c  (note: this differs from UMAT/Standard which uses 12,13,23)
c=====================================================================
      subroutine mat_from_vec(v, A)
      include 'vaba_param.inc'
      double precision v(6), A(3,3)
      A(1,1)=v(1)
      A(2,2)=v(2)
      A(3,3)=v(3)
      A(1,2)=v(4)
      A(2,1)=v(4)
      A(2,3)=v(5)
      A(3,2)=v(5)
      A(1,3)=v(6)
      A(3,1)=v(6)
      return
      end

c=====================================================================
c  HELPER: symmetric 3x3 tensor -> Voigt vector (6)
c  VUMAT order: [11, 22, 33, 12, 23, 13]
c=====================================================================
      subroutine vec_from_mat(A, v)
      include 'vaba_param.inc'
      double precision A(3,3), v(6)
      v(1)=A(1,1)
      v(2)=A(2,2)
      v(3)=A(3,3)
      v(4)=A(1,2)
      v(5)=A(2,3)
      v(6)=A(1,3)
      return
      end

c=====================================================================
c  HELPER: Zero a 3x3 tensor
c=====================================================================
      subroutine zero_tensor(T)
      include 'vaba_param.inc'
      double precision T(3,3)
      integer i,j
      do i=1,3
        do j=1,3
          T(i,j)=0.d0
        end do
      end do
      return
      end

