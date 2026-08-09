# Coupled Crystal Plasticity–Martensitic Transformation–Twinning VUMAT for TRIP Steels

**Abaqus/Explicit user material subroutine for transformation-induced plasticity, mechanical twinning, and crystallographic slip in metastable austenitic/TRIP steels**

---

## 1. Executive Summary

This repository contains a **VUMAT** user material subroutine for **Abaqus/Explicit** that implements a coupled constitutive model for **TRIP steels** and related metastable austenitic alloys. The model is designed to capture the simultaneous contribution of:

1. **Crystallographic slip** in FCC austenite using a rate-dependent crystal plasticity formulation.
2. **Martensitic transformation** from metastable retained austenite to martensite during plastic deformation.
3. **Mechanical twinning** on FCC twin systems.
4. **Twin–martensite interaction**, where twin formation can accelerate martensite nucleation.
5. **Transformation-induced volumetric expansion**, represented through an isotropic transformation eigenstrain.
6. **Numerical stabilization features** required for robust explicit finite element simulations.

The subroutine is written for a **homogenized material-point description**. Austenite is treated as the crystallographically active FCC phase carrying slip and twinning. Martensite is represented primarily by its **volume fraction**, its **volumetric transformation strain**, and its indirect contribution to the macroscopic hardening response. A full two-phase FCC/BCT crystal plasticity treatment is not implemented in the current version but is a natural extension.

The model is suitable for:

- TRIP steel forming simulations.
- Crashworthiness and energy-absorption studies.
- Polycrystalline RVE simulations when combined with grain orientations.
- Investigations of phase fraction evolution under mechanical loading.
- ICME workflows linking microstructure evolution to macroscopic mechanical response.

> **Important modeling note:**
> The theoretical framework behind this work can be formulated in a finite-strain multiplicative kinematic setting. The present VUMAT, however, uses an **incrementally additive, corotational/hypoelastic-style implementation** compatible with Abaqus/Explicit. It is robust and efficient for engineering simulations where elastic strains remain small compared with plastic, transformation, and twinning strains.

---

## 2. Physical Mechanisms Represented

| Physical mechanism | Representation in the VUMAT | Primary state variables |
|---|---|---|
| Elastic response | Isotropic Hookean elasticity | Stress tensor |
| Crystallographic slip | 12 FCC `{111}<110>` slip systems, rate-dependent power-law flow, Voce hardening | Slip shear strains, CRSS values |
| Martensitic transformation | Olson–Cohen-type shear-band nucleation + JMAK-type strain-driven kinetics | Martensite volume fraction $f_M$ |
| Transformation eigenstrain | Isotropic volumetric expansion proportional to $\Delta f_M$ | Accumulated transformation strain $\varepsilon_M$ |
| Mechanical twinning | 12 FCC `{111}<112>` twin systems, polar CRSS activation, TANH kinetics | Twin volume fraction $f_{\text{tw}}$ |
| Twin-induced plastic strain | Twin shear distributed over active twin systems | Twin plastic strain tensor |
| Twin–martensite coupling | Additional martensite nucleation proportional to twin fraction increment | Coupling parameter $\alpha_{TM}$ |
| Phase fraction constraint | $f_M + f_{\text{tw}} \leq 1$ enforced | $f_M$, $f_{\text{tw}}$ |
| Numerical stabilization | Slip increment cap, volume-fraction rate caps, anti-overshoot limiter | Internal algorithmic variables |

---

## 3. Constitutive Architecture

The model belongs to the class of **mechanism-based crystal plasticity models with internal phase evolution**. Instead of using a macroscopic yield surface such as von Mises, Hill48, or Barlat, the onset of plastic flow is governed locally on each crystallographic system by the resolved shear stress and the corresponding critical resolved shear stress.

The overall inelastic strain rate is conceptually decomposed as

$
\dot{\varepsilon}^{\text{inel}}
=
\dot{\varepsilon}^{\text{slip}}
+
\dot{\varepsilon}^{\text{twin}}
+
\dot{\varepsilon}^{\text{tr}}
$

where

- $\dot{\varepsilon}^{\text{slip}}$ is the plastic strain rate from dislocation glide,
- $\dot{\varepsilon}^{\text{twin}}$ is the pseudo-plastic strain rate associated with mechanical twinning,
- $\dot{\varepsilon}^{\text{tr}}$ is the transformation strain rate associated with martensite formation.

In the current implementation, the strain increment passed by Abaqus is treated incrementally:

$$
\Delta \varepsilon
=
\Delta \varepsilon^{e}
+
\Delta \varepsilon^{\text{slip}}
+
\Delta \varepsilon^{\text{twin}}
+
\Delta \varepsilon^{\text{tr,vol}}
$$

The stress update is then performed from the elastic part only:

$$
\Delta \sigma
=
\mathbf{D}^{e} :
\left(
\Delta \varepsilon
-
\Delta \varepsilon^{\text{slip}}
-
\Delta \varepsilon^{\text{twin}}
-
\Delta \varepsilon^{\text{tr,vol}}
\right)
$$

where $\mathbf{D}^{e}$ is the isotropic elastic stiffness tensor.

---

## 4. Kinematic Framework

### 4.1 Conceptual finite-strain picture

A rigorous finite-strain description of TRIP/TWIP behavior often uses a multiplicative decomposition of the deformation gradient:

$$
\mathbf{F}
=
\mathbf{F}^{e}
\mathbf{F}^{p}
\mathbf{F}^{\text{tr}}
\mathbf{F}^{\text{tw}}
$$

where

- $\mathbf{F}^{e}$ is the elastic lattice distortion,
- $\mathbf{F}^{p}$ is the plastic deformation due to crystallographic slip,
- $\mathbf{F}^{\text{tr}}$ is the transformation deformation due to martensite formation,
- $\mathbf{F}^{\text{tw}}$ is the twinning deformation.

This decomposition is useful because it separates physically distinct mechanisms:

| Component | Mechanism |
|---|---|
| $\mathbf{F}^{e}$ | Elastic stretching and lattice rotation |
| $\mathbf{F}^{p}$ | Dislocation slip |
| $\mathbf{F}^{\text{tr}}$ | Martensitic transformation strain |
| $\mathbf{F}^{\text{tw}}$ | Twinning shear and reorientation |

### 4.2 Implemented incremental form

The present VUMAT does not explicitly track all deformation gradients. Instead, it uses an **incremental strain decomposition** in the corotational frame supplied by Abaqus/Explicit:

$$
\Delta \varepsilon
=
\Delta \varepsilon^{e}
+
\Delta \varepsilon^{\text{cp}}
+
\Delta \varepsilon^{\text{tw}}
+
\Delta \varepsilon^{\text{tr,vol}}
$$

This approach is computationally efficient and robust for explicit simulations. It is particularly appropriate when:

- elastic strains are small,
- plastic strains are large,
- phase transformation strains are moderate,
- the main interest is macroscopic stress–strain response and phase fraction evolution.

For problems requiring exact finite-strain hyperelastic consistency, strong elastic anisotropy, or explicit lattice reorientation, a full multiplicative implementation would be preferable.

---

## 5. Elastic Response

### 5.1 Implemented model

The subroutine uses **isotropic linear elasticity**:

$$
\lambda = \frac{E \nu}{(1+\nu)(1-2\nu)}
$$

$$
\mu = \frac{E}{2(1+\nu)}
$$

$$
K = \lambda + \frac{2}{3}\mu
$$

The elastic stiffness matrix in Voigt notation is assembled as:

$$
\mathbf{D}^{e}
=
\begin{bmatrix}
\lambda + 2\mu & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda + 2\mu & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda + 2\mu & 0 & 0 & 0 \\
0 & 0 & 0 & 2\mu & 0 & 0 \\
0 & 0 & 0 & 0 & 2\mu & 0 \\
0 & 0 & 0 & 0 & 0 & 2\mu
\end{bmatrix}
$$

The stress increment is

$$
\Delta \sigma = \mathbf{D}^{e} : \Delta \varepsilon^{e}
$$

### 5.2 Elasticity model alternatives

| Elastic model | Advantages | Disadvantages | Suitability for present code |
|---|---|---|---|
| Isotropic Hooke | Simple, robust, only two parameters | Ignores crystal elastic anisotropy | Selected |
| Saint Venant–Kirchhoff | Finite-strain compatible | Poor for large strains, not ideal for metals | Not selected |
| Neo-Hookean hyperelastic | Thermodynamically consistent for finite elastic strains | Still isotropic, more complex | Possible extension |
| Cubic anisotropic elasticity | Captures FCC elastic anisotropy | Requires $C_{11}, C_{12}, C_{44}$ | Future extension |
| Anisotropic hyperelasticity | General texture-dependent elasticity | Many parameters, calibration burden | Not selected |

### 5.3 Why isotropic elasticity was selected

The current model focuses on the dominant sources of inelasticity in TRIP steels:

- dislocation slip,
- martensitic transformation,
- mechanical twinning,
- transformation/twinning-induced hardening.

For many engineering forming and impact simulations, the elastic strain is small compared with the plastic and transformation strains. Therefore, isotropic elasticity provides a good compromise between:

- numerical stability,
- computational efficiency,
- ease of calibration,
- robustness in explicit simulations.

The main plastic anisotropy is captured through the crystallographic slip and twin systems, not through elastic anisotropy.

---

## 6. Plastic Flow: Why Crystal Plasticity Instead of a Macroscopic Yield Criterion?

Several macroscopic yield criteria could be used to describe plastic yielding in steels:

| Yield/plasticity model | Strengths | Weaknesses |
|---|---|---|
| von Mises $J_2$ | Simple, robust, isotropic | Cannot capture crystallographic texture or orientation effects |
| Tresca | Simple shear-based criterion | Corners in yield surface, less accurate for FCC metals |
| Drucker–Prager | Includes hydrostatic pressure sensitivity | Mainly for soils, concrete, porous materials |
| Hill48 | Captures orthotropic sheet anisotropy | Empirical, no direct slip-system information |
| Barlat | Advanced sheet anisotropy | More parameters, higher complexity |
| Crystal plasticity | Directly models slip systems, texture, orientation, CRSS evolution | Higher computational cost |

### 6.1 Selected model: rate-dependent crystal plasticity

For TRIP steels, the physically dominant plastic mechanism in austenite is slip on FCC systems:

$$
\{111\}\langle 110 \rangle
$$

There are 12 such slip systems. Crystal plasticity is selected because it can represent:

- crystallographic orientation,
- texture evolution tendencies,
- grain-level anisotropy,
- slip-system-level hardening,
- coupling with twinning,
- coupling with phase transformation,
- localized deformation in polycrystalline RVEs.

In a crystal plasticity framework, the macroscopic yield surface is not imposed directly. Instead, yielding emerges from the collective activation of slip systems according to Schmid-type resolved shear stress conditions.

---

## 7. Crystallographic Slip Model

### 7.1 Slip systems

The model uses the standard FCC slip family:

$$
\{111\}\langle 110 \rangle
$$

with 12 systems:

| Plane family | Slip directions | Number of systems |
|---|---|---:|
| $(111)$ | $[0\bar{1}1]$, $[10\bar{1}]$, $[\bar{1}10]$ | 3 |
| $(1\bar{1}1)$ | $[011]$, $[101]$, $[\bar{1}\bar{1}0]$ | 3 |
| $(\bar{1}11)$ | $[01\bar{1}]$, $[101]$, $[\bar{1}10]$ | 3 |
| $(11\bar{1})$ | $[011]$, $[10\bar{1}]$, $[\bar{1}\bar{1}0]$ | 3 |

The slip systems are hard-coded in the local material coordinate system. For polycrystalline simulations, grain orientation must be supplied through Abaqus `*ORIENTATION`.

### 7.2 Schmid tensor

For each slip system $\alpha$, with slip direction $\mathbf{s}^{\alpha}$ and slip-plane normal $\mathbf{n}^{\alpha}$, the symmetric Schmid tensor is

$$
\mathbf{M}^{\alpha}
=
\frac{1}{2}
\left(
\mathbf{s}^{\alpha} \otimes \mathbf{n}^{\alpha}
+
\mathbf{n}^{\alpha} \otimes \mathbf{s}^{\alpha}
\right)
$$

The resolved shear stress is

$$
\tau^{\alpha}
=
\mathbf{M}^{\alpha} : \boldsymbol{\sigma}
$$

### 7.3 Rate-dependent power-law flow rule

The shear rate on slip system $\alpha$ is

$$
\dot{\gamma}^{\alpha}
=
\dot{\gamma}_{0}
\left|
\frac{\tau^{\alpha}}{g^{\alpha}}
\right|^{1/m}
\operatorname{sign}(\tau^{\alpha})
$$

where

- $\dot{\gamma}_{0}$ is the reference shear strain rate,
- $m$ is the rate-sensitivity exponent,
- $g^{\alpha}$ is the current critical resolved shear stress, CRSS.

The slip increment is

$$
\Delta \gamma^{\alpha}
=
\dot{\gamma}^{\alpha} \Delta t
$$

### 7.4 Plastic strain from slip

The plastic strain increment due to slip is

$$
\Delta \boldsymbol{\varepsilon}^{\text{cp}}
=
\sum_{\alpha=1}^{12}
\mathbf{M}^{\alpha}
\Delta \gamma^{\alpha}
$$

### 7.5 Voce hardening law

The CRSS of each slip system evolves according to a Voce-type saturation law:

$$
g^{\alpha}(\Gamma^{\alpha})
=
\tau_{0}
+
(\tau_{s} - \tau_{0})
\left[
1 -
\exp\left(
-\frac{h_{0}\Gamma^{\alpha}}{\tau_{s} - \tau_{0}}
\right)
\right]
$$

where

- $\tau_{0}$ is the initial CRSS,
- $\tau_{s}$ is the saturation CRSS,
- $h_{0}$ is the initial hardening rate,
- $\Gamma^{\alpha}$ is the accumulated absolute slip on system $\alpha$.

If $\tau_{s} \leq \tau_{0}$, the model falls back to linear hardening:

$$
g^{\alpha}
=
\tau_{0}
+
h_{0}\Gamma^{\alpha}
$$

### 7.6 Hardening model alternatives

| Hardening model | Advantages | Disadvantages | Reason for selection/rejection |
|---|---|---|---|
| Linear hardening | Simple | No saturation | Fallback only |
| Power-law hardening | Common for metals | Less natural saturation behavior | Not selected |
| Voce hardening | Captures saturation typical of FCC metals | Requires $\tau_s, h_0$ | Selected |
| Isotropic macroscopic hardening | Simple continuum model | No slip-system resolution | Not suitable for CP |
| Kinematic hardening | Captures Bauschinger effect | More parameters, not system-based here | Not implemented |
| Latent hardening matrix | Captures cross-system interaction | Requires latent hardening coefficients | Future extension |

### 7.7 Important slip-hardening simplification

The current implementation uses **self-hardening only**. The CRSS of each slip system evolves according to its own accumulated slip:

$$
\Gamma^{\alpha}
=
\sum |\Delta \gamma^{\alpha}|
$$

A more advanced model would include a latent hardening matrix:

$$
\dot{g}^{\alpha}
=
\sum_{\beta}
h^{\alpha\beta}
|\dot{\gamma}^{\beta}|
$$

Latent hardening is often important for texture evolution and multislip behavior in FCC metals. It is a recommended extension.

---

## 8. Martensitic Transformation

Martensitic transformation is one of the central mechanisms in TRIP steels. Metastable retained austenite transforms into martensite during plastic deformation, producing:

- additional hardening,
- delayed necking,
- increased strength,
- improved energy absorption,
- local volumetric expansion,
- transformation-induced plastic strain.

The transformation model in this VUMAT is built from three parts:

1. **Transformation kinetics**,
2. **Transformation eigenstrain**,
3. **Coupling with twinning and plastic strain**.

---

## 9. Thermodynamic Background of Martensitic Transformation

Although the present VUMAT uses a robust strain-driven implementation, the physical basis of martensitic transformation is thermodynamic.

### 9.1 Chemical driving force

The chemical free-energy difference between austenite and martensite can be written approximately as

$$
\Delta G_{\text{chem}}
=
\Delta H_{M}
-
T \Delta S_{M}
$$

where

- $\Delta H_{M}$ is the enthalpy change,
- $\Delta S_{M}$ is the entropy change,
- $T$ is absolute temperature.

Transformation becomes favorable when the total driving force exceeds a critical energy barrier.

### 9.2 Mechanical driving force

An applied stress state can assist transformation by performing mechanical work on the transformation strain:

$$
\Delta G_{\text{mech}}
=
\boldsymbol{\sigma} : \boldsymbol{\varepsilon}^{\text{tr}}
$$

Thus the total driving force may be written conceptually as

$$
\Delta G_{M}
=
\Delta G_{\text{chem}}
+
\Delta G_{\text{mech}}
$$

or, in a more explicit form,

$$
\Delta G_{M}
=
\Delta H_{M}
-
T \Delta S_{M}
+
\boldsymbol{\sigma} : \boldsymbol{\varepsilon}^{\text{tr}}
$$

### 9.3 Critical barrier and Olson–Cohen-type condition

In Olson–Cohen-type physical nucleation models, transformation starts when the driving force exceeds a critical barrier:

$$
\Delta G_{M} \geq G^{\text{crit}}_{M}
$$

A common simplified form is

$$
G^{\text{crit}}_{M}
=
G^{0}_{M}
-
c \sigma^{2}
$$

where

- $G^{0}_{M}$ is the stress-free barrier,
- $c$ is a material constant,
- $\sigma$ represents the applied stress intensity.

This expression captures the idea that mechanical loading lowers the nucleation barrier.

### 9.4 Stress-assisted versus strain-induced transformation

Two regimes are commonly distinguished:

| Regime | Description |
|---|---|
| Stress-assisted transformation | Transformation is assisted by applied stress before large plastic strain; often important near or below $M_s$ |
| Strain-induced transformation | Transformation is driven by plastic deformation, shear bands, and defect generation; dominant in many TRIP steels |

The present VUMAT primarily captures the **strain-induced transformation regime** through accumulated plastic strain and shear-band-like nucleation variables.

### 9.5 Stress triaxiality and Lode angle

Advanced transformation models may include the effect of:

- stress triaxiality,
- Lode angle,
- hydrostatic pressure,
- shear-dominant versus tension-dominant loading.

Stress triaxiality is often defined as

$$
\eta
=
\frac{\sigma_m}{\bar{\sigma}}
$$

where

$$
\sigma_m = \frac{1}{3}\operatorname{tr}(\boldsymbol{\sigma})
$$

and $\bar{\sigma}$ is the von Mises equivalent stress.

High positive triaxiality generally promotes volumetric expansion associated with martensitic transformation, while shear-dominant states may alter the transformation rate and variant selection.

The current VUMAT does not explicitly compute triaxiality or Lode angle inside the kinetic law. However, these effects can be introduced in future versions by modifying the transformation kinetics parameters as functions of the local stress state.

### 9.6 Temperature and adiabatic heating

Martensitic transformation is temperature dependent. Increasing temperature generally stabilizes austenite and reduces the transformation rate. At high strain rates, adiabatic heating may raise the local temperature and suppress transformation.

The current implementation stores temperature as a state variable but keeps it fixed by default. Temperature dependence can be introduced externally by making parameters such as:

- `K_J`,
- `tau_twin`,
- `tau0`,
- `taus`,
- `beta_OC`,

functions of temperature.

---

## 10. Martensitic Transformation Kinetics

### 10.1 General kinetic behavior

The evolution of martensite volume fraction versus plastic strain in TRIP steels is often sigmoidal:

1. an initial incubation stage,
2. a rapid transformation stage,
3. a saturation stage as austenite is consumed.

A robust kinetic model should capture these three stages.

The present model combines two kinetic descriptions:

1. **Olson–Cohen-type shear-band nucleation**,
2. **JMAK-type transformation kinetics**.

They are combined as independent transformation mechanisms.

---

## 11. Olson–Cohen-Type Shear-Band Nucleation Model

The Olson–Cohen model is one of the most physically meaningful models for strain-induced martensitic transformation. It assumes that plastic deformation generates shear bands, and intersections of shear bands act as martensite nucleation sites.

### 11.1 Shear-band density evolution

The VUMAT uses a normalized shear-band density variable $N_{\text{sb}}$. Its increment is

$$
\Delta N_{\text{sb}}
=
\alpha_{\text{OC}}
(1 - N_{\text{sb}})
\frac{\Delta \bar{\varepsilon}^{p}}{\varepsilon_{0,\text{OC}}}
$$

where

- $\alpha_{\text{OC}}$ controls the shear-band generation rate,
- $\varepsilon_{0,\text{OC}}$ is a reference plastic strain,
- $\Delta \bar{\varepsilon}^{p}$ is the equivalent plastic strain increment.

The updated shear-band density is capped between 0 and 1:

$$
0 \leq N_{\text{sb}} \leq 1
$$

### 11.2 Martensite fraction from Olson–Cohen nucleation

The transformed fraction associated with shear-band nucleation is

$$
f_{M}^{\text{OC}}
=
1 -
\exp
\left[
-\beta_{\text{OC}}
N_{\text{sb}}^{n_{\text{OC}}}
\right]
$$

where

- $\beta_{\text{OC}}$ controls nucleation efficiency,
- $n_{\text{OC}}$ controls the shape of the transformation curve.

### 11.3 Advantages of the Olson–Cohen model

- Physically based on shear-band nucleation.
- Captures strain-induced transformation.
- Naturally produces sigmoidal behavior.
- Compatible with crystal plasticity.
- Parameters have microstructural interpretation.

### 11.4 Limitations

- Does not explicitly include stress triaxiality or Lode angle.
- Does not explicitly include temperature unless parameters are modified.
- Does not distinguish individual martensite variants.
- Assumes austenite transformation is controlled primarily by plastic strain and shear-band density.

---

## 12. JMAK-Type Transformation Kinetics

The Johnson–Mehl–Avrami–Kolmogorov model is a classical nucleation-and-growth model. In strain-driven form, it may be written as

$$
f_{M}^{\text{JMAK}}
=
1 -
\exp
\left[
-K_{J}
\left(
\bar{\varepsilon}^{p}
\right)^{n_{J}}
\right]
$$

where

- $K_{J}$ is a temperature-dependent rate constant,
- $n_{J}$ is the Avrami exponent,
- $\bar{\varepsilon}^{p}$ is the accumulated equivalent plastic strain.

### 12.1 Advantages of JMAK

- Captures sigmoidal transformation behavior.
- Simple and robust.
- Useful for fitting experimental transformation curves.
- Can absorb temperature dependence through $K_J$.

### 12.2 Limitations

- Less mechanistic than Olson–Cohen.
- Does not explicitly represent shear-band nucleation.
- Does not directly include stress-state effects.
- Requires calibration against phase fraction data.

---

## 13. Combined Olson–Cohen and JMAK Kinetics

The VUMAT combines the two transformation descriptions as independent mechanisms. If two independent transformation mechanisms produce fractions $f_{M}^{\text{OC}}$ and $f_{M}^{\text{JMAK}}$, the combined fraction is

$$
f_{M}^{\text{base}}
=
1 -
\left(
1 - f_{M}^{\text{OC}}
\right)
\left(
1 - f_{M}^{\text{JMAK}}
\right)
$$

This expression means that the untransformed austenite fraction is the product of the untransformed fractions from each mechanism:

$$
1 - f_{M}^{\text{base}}
=
\left(
1 - f_{M}^{\text{OC}}
\right)
\left(
1 - f_{M}^{\text{JMAK}}
\right)
$$

### 13.1 Why combine OC and JMAK?

| Model | Captures |
|---|---|
| Olson–Cohen | Shear-band nucleation, strain-induced transformation |
| JMAK | Overall sigmoidal nucleation/growth behavior |
| Combined model | More flexible representation of experimental data |

The combined model allows the user to represent materials where transformation is controlled by both defect generation and overall nucleation/growth statistics.

---

## 14. Twin–Martensite Coupling

Mechanical twins can act as additional nucleation sites for martensite. To represent this, the model adds a coupling term:

$$
f_{M}^{\text{base}}
\leftarrow
f_{M}^{\text{base}}
+
\alpha_{TM}
\Delta f_{\text{tw}}
\left(
1 - f_{M}^{\text{base}}
\right)
$$

where

- $\alpha_{TM}$ is the twin-to-martensite coupling coefficient,
- $\Delta f_{\text{tw}}$ is the twin volume fraction increment,
- $1 - f_{M}^{\text{base}}$ is the remaining transformable austenite.

If $\alpha_{TM}=0$, the coupling is disabled.

### 14.1 Physical meaning

This term captures the idea that twin boundaries can:

- act as martensite nucleation sites,
- increase local defect density,
- assist transformation in metastable austenite,
- alter local transformation kinetics.

### 14.2 Monotonicity constraint

Martensite transformation is treated as irreversible:

$$
f_{M}^{\text{new}} \geq f_{M}^{\text{old}}
$$

This prevents nonphysical reversal of martensite to austenite during unloading.

---

## 15. Transformation Eigenstrain

Martensitic transformation in steels is accompanied by a volume expansion. The present model represents this using an isotropic volumetric eigenstrain.

Let

$$
\Delta V/V
$$

be the volumetric transformation strain associated with complete transformation. In the code this parameter is `dVoV`. Typical values for austenite-to-martensite transformation are approximately

$$
0.02 \leq \Delta V/V \leq 0.04
$$

For a martensite fraction increment $\Delta f_M$, the volumetric eigenstrain increment per spatial direction is

$$
\Delta \varepsilon^{\text{tr}}_{ii}
=
\frac{\Delta f_M \Delta V/V}{3}
$$

The corresponding stress correction for constrained expansion is

$$
\Delta \sigma_{ii}
=
-K \Delta f_M \Delta V/V
$$

for the three normal components.

### 15.1 Transformation effect model alternatives

| Transformation strain model | Description | Advantages | Disadvantages |
|---|---|---|---|
| Unity model | No transformation strain | Simplest | Physically inadequate for TRIP steels |
| Approximate isotropic volumetric model | Small-strain volumetric expansion | Simple | Less exact for finite strain |
| Exact volumetric model | Uses transformation Jacobian $J^{\text{tr}}$ | Thermodynamically cleaner | Still isotropic |
| Directional shear model | Includes variant-specific shear | Captures variant selection | Requires variant data and calibration |
| Combined volume-shear model | Includes dilatation and shear | Most complete | Complex, many parameters |
| Statistical/microstructural model | Spatially heterogeneous transformation | Captures local effects | High computational cost |

### 15.2 Why the volumetric eigenstrain model was selected

The current model uses an isotropic volumetric transformation strain because:

1. The dominant macroscopic effect of austenite-to-martensite transformation is volume expansion.
2. The model remains robust for explicit simulations.
3. It requires only one main parameter, `dVoV`.
4. It avoids the need to track individual martensite variants.
5. It can be calibrated from dilatometry or crystallographic data.

A directional shear or combined volume-shear transformation model is a recommended extension when variant selection, transformation texture, or shear-assisted transformation is important.

---

## 16. Mechanical Twinning

Mechanical twinning is another important deformation mechanism in low-SFE FCC alloys, TWIP steels, and some TRIP-assisted steels. Twinning contributes to:

- plastic strain,
- work hardening,
- texture evolution,
- dynamic Hall–Petch-type strengthening,
- martensite nucleation at twin boundaries.

The present model represents twinning using 12 FCC twin systems:

$$
\{111\}\langle 112 \rangle
$$

---

## 17. Twinning Kinematic Model

### 17.1 Twin Schmid tensors

For each twin system $t$, a symmetric Schmid-like tensor is constructed:

$$
\mathbf{M}^{t}
=
\frac{1}{2}
\left(
\mathbf{s}^{t} \otimes \mathbf{n}^{t}
+
\mathbf{n}^{t} \otimes \mathbf{s}^{t}
\right)
$$

where

- $\mathbf{n}^{t}$ is the twin-plane normal,
- $\mathbf{s}^{t}$ is the twin shear direction.

The resolved twin shear stress is

$$
\tau^{t}
=
\mathbf{M}^{t} : \boldsymbol{\sigma}
$$

### 17.2 Polar CRSS activation

Twinning is inherently directional. Unlike ordinary slip, which can often be treated as approximately bidirectional, twinning occurs in a specific crystallographic sense.

The model therefore uses a **polar CRSS criterion**:

$$
\tau^{t} > \tau_{c}^{\text{tw}}
$$

and also requires

$$
\tau^{t} > 0
$$

Thus a twin system activates only when the resolved shear stress is both positive and larger than the critical twin shear stress.

This is more physical than a non-polar absolute-value criterion.

### 17.3 Twinning condition model alternatives

| Twinning criterion | Advantages | Disadvantages |
|---|---|---|
| Absolute CRSS | Simple | Does not distinguish twin sense |
| Polar CRSS | Captures directional twinning | Requires sign convention care |
| Energy-based SFE criterion | Thermodynamically richer | Requires stacking fault energy data |
| Full nucleation theory | Most physical | Complex and data-intensive |

The selected model is the **polar CRSS criterion** because it is simple, robust, and physically adequate for engineering-scale simulations.

### 17.4 Twin systems

| Plane | Shear directions ($\pm \langle 112 \rangle/\sqrt{6}$) | Systems |
|---|---|---|
| $(1\ 1\ 1)$ | $[1\ 1\ \bar{2}]$, $[1\ \bar{2}\ 1]$, $[\bar{2}\ 1\ 1]$ | 1–3 |
| $(1\ \bar{1}\ 1)$ | $[2\ 1\ \bar{1}]$, $[1\ 2\ 1]$, $[\bar{1}\ 1\ 2]$ | 4–6 |
| $(\bar{1}\ 1\ 1)$ | $[1\ 2\ \bar{1}]$, $[2\ 1\ 1]$, $[1\ \bar{1}\ 2]$ | 7–9 |
| $(1\ 1\ \bar{1})$ | $[1\ 1\ 2]$, $[2\ \bar{1}\ 1]$, $[\bar{1}\ 2\ 1]$ | 10–12 |

$\mathbf{n}\cdot\mathbf{s}=0$ holds for every pair; the sign of $\mathbf{s}$ sets the twinning sense (verify against your crystallography).

---

## 18. Twinning Kinetics

### 18.1 TANH kinetic law

The model uses a smooth hyperbolic tangent law to determine the target twin volume fraction:

$$
f_{\text{tw}}^{\text{target}}
=
A_{\text{tw}}
\frac{1}{2}
\left[
1 +
\tanh
\left(
\frac{
\tau_{\max}^{t}
-
\tau_{c}^{\text{tw}}
}{
B_{\text{tw}}
}
\right)
\right]
$$

where

- $A_{\text{tw}}$ is the maximum twin volume fraction,
- $B_{\text{tw}}$ controls the transition width,
- $\tau_{\max}^{t}$ is the maximum resolved twin shear stress among all twin systems,
- $\tau_{c}^{\text{tw}}$ is the critical twin shear stress.

### 18.2 Irreversibility

Twin growth is treated as irreversible:

$$
\Delta f_{\text{tw}}
=
\max
\left(
0,
f_{\text{tw}}^{\text{target}} - f_{\text{tw}}^{\text{old}}
\right)
$$

This prevents detwinning during unloading in the current implementation.

### 18.3 Rate cap for numerical stability

The twin fraction increment is limited by a maximum value per increment:

$$
\Delta f_{\text{tw}} \leq \Delta f_{\text{tw}}^{\max}
$$

In the current code, the default cap is

$$
\Delta f_{\text{tw}}^{\max} = 0.02
$$

This improves explicit stability and prevents abrupt strain spikes.

### 18.4 Twin plastic strain

The twin plastic strain increment is distributed over active twin systems according to their relative driving stresses.

Define the active-system weight:

$$
w^{t}
=
\max
\left(
0,
\tau^{t} - \tau_{c}^{\text{tw}}
\right)
$$

The total weight is

$$
W =
\sum_{t} w^{t}
$$

If $W > 0$, the twin plastic strain increment is

$$
\Delta \boldsymbol{\varepsilon}^{\text{tw}}
=
\sum_{t}
\frac{w^{t}}{W}
\Delta f_{\text{tw}}^{\text{eff}}
\gamma_{\text{tw}}
\mathbf{M}^{t}
$$

where

- $\gamma_{\text{tw}}$ is the characteristic twin shear,
- for FCC twinning, $\gamma_{\text{tw}} \approx 0.707$.

### 18.5 Twinning kinetics alternatives

| Kinetic model | Advantages | Disadvantages | Reason for selection/rejection |
|---|---|---|---|
| Linear kinetic law | Simple | No saturation, unrealistic at large strain | Not selected |
| Logistic kinetic law | Captures saturation | Requires rate parameters | Possible alternative |
| PTR, predominant twin reorientation | Captures texture reorientation | Numerically discontinuous, complex | Future extension |
| TANH kinetic law | Smooth, stable, continuous, captures gradual activation and saturation | Phenomenological | Selected |
| Fully coupled twin-slip crystal plasticity | Most physical | High cost and complexity | Future extension |

### 18.6 Why TANH kinetics were selected

The hyperbolic tangent law was selected because it provides:

- smooth activation,
- continuous derivatives,
- stable explicit integration,
- natural saturation,
- easy calibration,
- no abrupt jumps in volume fraction.

This is particularly important in Abaqus/Explicit, where sudden internal strain increments can destabilize the solution.

---

## 19. Phase Fraction Constraint

Because both martensite and twins consume austenite, the model enforces

$$
f_M + f_{\text{tw}} \leq 1
$$

The implementation proceeds as follows:

1. A preliminary twin increment is computed and capped by the available austenite.
2. Martensite transformation is computed, including twin-assisted coupling.
3. The twin fraction is re-checked against the updated martensite fraction.
4. If necessary, the twin fraction is reduced so that the total transformed fraction does not exceed unity.

This prevents nonphysical phase fractions.

---

## 20. Complete Incremental Algorithm

For each material point and each increment, the VUMAT performs the following steps.

### Step 0: Initialization

If `stepTime = 0`, Abaqus is performing an initial check. The subroutine performs a purely elastic update and does not activate plasticity, twinning, or transformation.

If `totalTime = 0`, state variables are initialized:

- slip strains = 0,
- CRSS = `tau0`,
- martensite fraction = 0,
- twin fraction = 0,
- equivalent plastic strain = 0,
- shear-band density = 0,
- transformation strain = 0,
- temperature = 298 K by default.

### Step 1: Elastic predictor

A trial stress is computed:

$$
\boldsymbol{\sigma}^{\text{trial}}
=
\boldsymbol{\sigma}^{\text{old}}
+
\mathbf{D}^{e} : \Delta \boldsymbol{\varepsilon}
$$

### Step 2: Slip-system resolved shear stresses

For each of the 12 slip systems:

$$
\tau^{\alpha}
=
\mathbf{M}^{\alpha} : \boldsymbol{\sigma}^{\text{trial}}
$$

### Step 3: Rate-dependent slip increments

$$
\Delta \gamma^{\alpha}
=
\dot{\gamma}_{0}
\left|
\frac{\tau^{\alpha}}{g^{\alpha}}
\right|^{1/m}
\operatorname{sign}(\tau^{\alpha})
\Delta t
$$

The slip increments are capped for stability.

### Step 4: Assemble slip plastic strain

$$
\Delta \boldsymbol{\varepsilon}^{\text{cp}}
=
\sum_{\alpha=1}^{12}
\mathbf{M}^{\alpha}
\Delta \gamma^{\alpha}
$$

### Step 5: Anti-overshoot limiter

The deviatoric part of the slip strain increment is compared with the deviatoric part of the imposed strain increment. If the slip strain increment is larger, it is scaled back:

$$
\Delta \boldsymbol{\varepsilon}^{\text{cp}}
\leftarrow
s
\Delta \boldsymbol{\varepsilon}^{\text{cp}}
$$

with

$$
s =
\frac{
\left\|
\operatorname{dev}
\Delta \boldsymbol{\varepsilon}
\right\|
}{
\left\|
\operatorname{dev}
\Delta \boldsymbol{\varepsilon}^{\text{cp}}
\right\|
}
$$

This prevents explicit instability caused by excessive plastic correction in a single increment.

### Step 6: Update Voce hardening

For each slip system:

$$
\Gamma^{\alpha}_{\text{new}}
=
\Gamma^{\alpha}_{\text{old}}
+
|\Delta \gamma^{\alpha}|
$$

Then $g^{\alpha}$ is updated using the Voce law.

### Step 7: Twinning resolved shear stresses

For each twin system:

$$
\tau^{t}
=
\mathbf{M}^{t} : \boldsymbol{\sigma}^{\text{trial}}
$$

### Step 8: TANH twin target fraction

$$
f_{\text{tw}}^{\text{target}}
=
A_{\text{tw}}
\frac{1}{2}
\left[
1 +
\tanh
\left(
\frac{
\tau_{\max}^{t}
-
\tau_{c}^{\text{tw}}
}{
B_{\text{tw}}
}
\right)
\right]
$$

### Step 9: Twin fraction update

The twin fraction increment is made irreversible and capped.

### Step 10: Martensite transformation

The equivalent plastic strain increment is computed from the slip strain increment:

$$
\Delta \bar{\varepsilon}^{p}
=
\sqrt{
\frac{2}{3}
\Delta \boldsymbol{\varepsilon}^{\text{cp}}
:
\Delta \boldsymbol{\varepsilon}^{\text{cp}}
}
$$

Then the Olson–Cohen shear-band density is updated:

$$
N_{\text{sb,new}}
=
N_{\text{sb,old}}
+
\alpha_{\text{OC}}
(1 - N_{\text{sb,old}})
\frac{
\Delta \bar{\varepsilon}^{p}
}{
\varepsilon_{0,\text{OC}}
}
$$

The Olson–Cohen martensite fraction is

$$
f_{M}^{\text{OC}}
=
1 -
\exp
\left[
-\beta_{\text{OC}}
N_{\text{sb,new}}^{n_{\text{OC}}}
\right]
$$

The JMAK fraction is

$$
f_{M}^{\text{JMAK}}
=
1 -
\exp
\left[
-K_{J}
\left(
\bar{\varepsilon}^{p}_{\text{new}}
\right)^{n_{J}}
\right]
$$

The combined base fraction is

$$
f_{M}^{\text{base}}
=
1 -
\left(
1 - f_{M}^{\text{OC}}
\right)
\left(
1 - f_{M}^{\text{JMAK}}
\right)
$$

Twin coupling is then added:

$$
f_{M}^{\text{base}}
\leftarrow
f_{M}^{\text{base}}
+
\alpha_{TM}
\Delta f_{\text{tw}}
\left(
1 - f_{M}^{\text{base}}
\right)
$$

The martensite fraction is enforced to be monotonic and bounded:

$$
f_{M}^{\text{old}}
\leq
f_{M}^{\text{new}}
\leq
1
$$

The increment is also capped for explicit stability.

### Step 11: Re-check twin fraction

The twin fraction is reduced if necessary to satisfy

$$
f_{M}^{\text{new}} + f_{\text{tw}}^{\text{new}} \leq 1
$$

### Step 12: Assemble total inelastic strain

$$
\Delta \boldsymbol{\varepsilon}^{\text{inel}}
=
\Delta \boldsymbol{\varepsilon}^{\text{cp}}
+
\Delta \boldsymbol{\varepsilon}^{\text{tw}}
$$

### Step 13: Stress update

The elastic strain increment is

$$
\Delta \boldsymbol{\varepsilon}^{e}
=
\Delta \boldsymbol{\varepsilon}
-
\Delta \boldsymbol{\varepsilon}^{\text{inel}}
$$

The stress increment is

$$
\Delta \boldsymbol{\sigma}
=
\mathbf{D}^{e}
:
\Delta \boldsymbol{\varepsilon}^{e}
$$

Then the martensite volumetric eigenstrain correction is applied:

$$
\Delta \sigma_{ii}
\leftarrow
\Delta \sigma_{ii}
-
K \Delta f_M \Delta V/V
$$

### Step 14: State variable update

All internal variables are written to `stateNew`.

---

## 21. Material Property Array

The subroutine expects:

    nprops = 20

| Index | Symbol | Meaning | Typical unit |
|---:|---|---|---|
| 1 | `E` | Young's modulus | stress unit, e.g. MPa |
| 2 | `nu` | Poisson's ratio | — |
| 3 | `gam0` | Reference slip rate | 1/s |
| 4 | `mexp` | Rate-sensitivity exponent | — |
| 5 | `tau0` | Initial CRSS | stress unit |
| 6 | `taus` | Saturation CRSS | stress unit |
| 7 | `h0` | Voce hardening modulus | stress unit |
| 8 | `alpha_OC` | Olson–Cohen shear-band generation coefficient | — |
| 9 | `beta_OC` | Olson–Cohen nucleation coefficient | — |
| 10 | `n_OC` | Olson–Cohen exponent | — |
| 11 | `eps0_OC` | Reference plastic strain for Olson–Cohen | — |
| 12 | `K_J` | JMAK rate constant | — |
| 13 | `n_J` | JMAK/Avrami exponent | — |
| 14 | `T_ref` | Reference temperature | K |
| 15 | `dVoV` | Volumetric transformation strain $\Delta V/V$ | — |
| 16 | `tau_twin` | Critical resolved shear stress for twinning | stress unit |
| 17 | `gam_twin` | Characteristic twin shear | — |
| 18 | `A_tw` | Maximum twin volume fraction | — |
| 19 | `B_tw` | TANH transition width | stress unit |
| 20 | `alpha_TM` | Twin-to-martensite coupling coefficient | — |

---

## 22. State Variables

The subroutine requires:

    nstatev = 30

| State variable | Index | Meaning |
|---:|---:|---|
| 1–12 | `sv(1:12)` | Accumulated slip shear strain on slip systems 1–12 |
| 13–24 | `sv(13:24)` | Current CRSS on slip systems 1–12 |
| 25 | `sv(25)` | Martensite volume fraction $f_M$ |
| 26 | `sv(26)` | Twin volume fraction $f_{\text{tw}}$ |
| 27 | `sv(27)` | Equivalent plastic strain from slip |
| 28 | `sv(28)` | Normalized shear-band density $N_{\text{sb}}$ |
| 29 | `sv(29)` | Accumulated transformation strain |
| 30 | `sv(30)` | Temperature |

These can be requested as Abaqus field outputs using `SDV1` through `SDV30`.

Recommended outputs:

| Output | SDV | Use |
|---|---:|---|
| Martensite fraction | 25 | Phase evolution |
| Twin fraction | 26 | TWIP evolution |
| Equivalent plastic strain | 27 | Plastic localization |
| Shear-band density | 28 | Transformation nucleation |
| Transformation strain | 29 | Volumetric transformation effect |
| Slip-system CRSS | 13–24 | Hardening evolution |

---

## 23. Required Abaqus Input Deck

A minimal material definition should contain:

    *MATERIAL, NAME=TRIP_STEEL_CP
    *DENSITY
    <density>
    *DEPVAR
    30
    *USER MATERIAL, CONSTANTS=20
    <E>, <nu>, <gam0>, <mexp>,
    <tau0>, <taus>, <h0>,
    <alpha_OC>, <beta_OC>, <n_OC>, <eps0_OC>,
    <K_J>, <n_J>, <T_ref>, <dVoV>,
    <tau_twin>, <gam_twin>, <A_tw>, <B_tw>, <alpha_TM>

Example placeholder:

    *MATERIAL, NAME=TRIP_STEEL_CP
    *DENSITY
    7.85e-9
    *DEPVAR
    30
    *USER MATERIAL, CONSTANTS=20
    200000., 0.3, 1.0, 0.05,
    120., 300., 1000.,
    10., 5., 4.5, 0.05,
    0.8, 2.0, 298., 0.03,
    180., 0.707, 0.25, 50., 0.5

The above values are placeholders only and must be calibrated for the specific steel grade.

---

## 24. Crystallographic Orientation

The slip and twin systems are defined in the local material coordinate system. Therefore, for polycrystalline simulations, each grain or grain group must be assigned an orientation.

Example:

    *ORIENTATION, NAME=GRAIN_1, SYSTEM=COORDINATES
    1., 0., 0., 0., 1., 0., 0., 0., 1.

For Euler-angle-based orientation assignment, use the appropriate Abaqus orientation definition or preprocessing tool.

For a Voronoi RVE:

- create one element set per grain,
- assign one orientation per grain,
- assign the same TRIP VUMAT material to all grains,
- let Abaqus rotate the material frame grain by grain.

The VUMAT itself does not store or rotate grain orientations internally.

---

## 25. Voigt Component Ordering

Abaqus/Explicit VUMAT uses tensorial shear strain components, not engineering shear strains.

The helper routines in this subroutine assume the following order:

    [11, 22, 33, 12, 23, 13]

This is different from some Abaqus/Standard UMAT conventions. Users modifying the tensor conversion routines should be careful.

---

## 26. Numerical Stability Features

Explicit simulations can become unstable if internal variables evolve too rapidly. This VUMAT includes several protective features.

### 26.1 Power-law cap

The power-law term is capped:

$$
\left|
\frac{\tau}{g}
\right|^{1/m}
\leq
10^{8}
$$

This prevents overflow for very high stress ratios.

### 26.2 Slip increment cap

Each slip increment is limited:

$$
|\Delta \gamma^{\alpha}| \leq 0.01
$$

This prevents excessively large slip in one increment.

### 26.3 Anti-overshoot limiter

The total deviatoric slip strain increment cannot exceed the imposed deviatoric strain increment. This avoids stress reversal due to over-correction.

### 26.4 Twin fraction rate cap

$$
\Delta f_{\text{tw}} \leq 0.02
$$

### 26.5 Martensite fraction rate cap

$$
\Delta f_{M} \leq 0.02
$$

### 26.6 Argument clamping

Exponential and hyperbolic tangent arguments are clamped to safe ranges to prevent overflow or numerical singularities.

### 26.7 Monotonic phase evolution

Martensite and twin fractions are prevented from decreasing during unloading in the current implementation.

---

## 27. Strengths of the Implemented Model

### 27.1 Mechanism-based description

The model directly represents the main deformation and transformation mechanisms in TRIP/TWIP steels:

- slip,
- twinning,
- martensitic transformation,
- transformation dilatation,
- twin-assisted transformation.

### 27.2 Crystallographic resolution

Because slip and twinning are defined on crystallographic systems, the model can represent orientation-dependent behavior and can be used in polycrystalline RVE simulations.

### 27.3 Robust explicit implementation

The subroutine includes several stabilization features that are essential for Abaqus/Explicit:

- slip caps,
- volume-fraction caps,
- anti-overshoot limiter,
- smooth TANH kinetics,
- bounded phase fractions.

### 27.4 Flexible transformation kinetics

The combination of Olson–Cohen and JMAK kinetics gives flexibility in fitting experimental phase fraction evolution.

### 27.5 Twin–martensite interaction

The coupling term allows the model to capture the experimentally observed effect that twins may promote martensite nucleation.

### 27.6 Practical calibratability

The model uses a manageable number of physically meaningful parameters. It can be calibrated using:

- uniaxial tension,
- XRD phase fraction measurements,
- EBSD texture data,
- dilatometry,
- digital image correlation,
- strain-rate sensitivity tests.

---

## 28. Limitations and Modeling Simplifications

The model is intentionally designed as a robust engineering model. Users should understand its simplifications.

### 28.1 Isotropic elasticity

Elastic anisotropy of FCC austenite is not included. This is acceptable for many engineering simulations but may affect texture-sensitive elastic response.

### 28.2 Self-hardening only

Latent hardening between slip systems is not currently included. This can affect multislip hardening and texture evolution.

### 28.3 Martensite is not a separate crystal plasticity phase

Martensite is represented by volume fraction and eigenstrain, not by its own BCT/BCC slip systems. Therefore, the model does not directly capture martensite plasticity or two-phase load partitioning at the crystal level.

### 28.4 No explicit variant selection

The current transformation eigenstrain is isotropic volumetric expansion. Martensite variant selection and transformation texture are not modeled.

### 28.5 No explicit stress triaxiality or Lode angle dependence

The current kinetic laws do not explicitly include stress triaxiality or Lode angle. These can be added in future versions.

### 28.6 Temperature is fixed by default

The model stores temperature but does not update it from plastic dissipation. Adiabatic heating is not automatically included.

### 28.7 Twinning is irreversible and unidirectional

Detwinning, twin reorientation, and predominant twin reorientation are not included.

### 28.8 Energy output is approximate

The internal energy and inelastic energy outputs are not fully thermodynamically consistent. They should not be used for precise energy balance studies unless the energy update is improved.

### 28.9 Forward explicit update

The constitutive update is not a fully implicit return-mapping algorithm. It is designed for Abaqus/Explicit with small stable increments.

---

## 29. Comparison with Macroscopic Constitutive Models

| Feature | Macroscopic $J_2$/Hill/Barlat model | Present crystal plasticity TRIP model |
|---|---|---|
| Crystallographic slip | No | Yes |
| Texture sensitivity | Limited or empirical | Yes |
| Slip-system hardening | No | Yes |
| Twin systems | No | Yes |
| Phase fraction evolution | Usually no | Yes |
| Transformation eigenstrain | Usually no | Yes |
| Twin–martensite coupling | No | Yes |
| Computational cost | Low | Moderate to high |
| Calibration complexity | Lower | Higher |
| Physical resolution | Low | High |

---

## 30. Recommended Calibration Strategy

### 30.1 Elastic parameters

Obtain from literature or experiments:

- $E$,
- $\nu$.

For steel, typical values are approximately:

$$
E \approx 190\text{--}210 \ \text{GPa}
$$

$$
\nu \approx 0.28\text{--}0.30
$$

Use consistent Abaqus units.

### 30.2 Slip parameters

Calibrate using single-crystal or polycrystalline stress–strain data:

- `tau0` from initial yield,
- `taus` from saturation stress,
- `h0` from initial hardening slope,
- `gam0` and `mexp` from strain-rate sensitivity.

For a rough initial estimate:

$$
\sigma_{\text{yield}}
\approx
\frac{\tau_0}{M_{\text{Schmid,max}}}
$$

where $M_{\text{Schmid,max}}$ is the maximum Schmid factor for the loading direction.

### 30.3 Transformation parameters

Use experimental martensite fraction versus plastic strain data.

Fit:

- `alpha_OC`,
- `beta_OC`,
- `n_OC`,
- `eps0_OC`,
- `K_J`,
- `n_J`.

The Olson–Cohen parameters mainly control the early transformation and shear-band nucleation behavior. The JMAK parameters control the overall sigmoidal evolution.

### 30.4 Volumetric transformation strain

Use crystallographic data or dilatometry.

Typical values:

$$
\Delta V/V \approx 0.02\text{--}0.04
$$

### 30.5 Twinning parameters

Calibrate using twin fraction measurements or mechanical response:

- `tau_twin` from twinning onset,
- `gam_twin` from crystallographic twin shear, usually about 0.707 for FCC,
- `A_tw` from maximum observed twin fraction,
- `B_tw` from the sharpness of twin activation.

### 30.6 Twin–martensite coupling

Set `alpha_TM = 0` initially. Then increase it if experimental evidence shows that twinning accelerates martensite formation.

---

## 31. Suggested Validation Cases

### 31.1 Single-element uniaxial tension

Check:

- elastic slope,
- initial yield,
- hardening shape,
- stability.

### 31.2 Single-crystal orientation test

Apply tension along known crystallographic directions and compare yield stress against Schmid-factor predictions.

### 31.3 Polycrystalline RVE

Assign random or measured orientations and examine:

- macroscopic stress–strain response,
- phase fraction evolution,
- strain localization,
- grain-level heterogeneity.

### 31.4 Phase fraction evolution

Compare predicted $f_M$ versus strain against:

- XRD,
- EBSD,
- magnetic measurements,
- neutron diffraction.

### 31.5 Twin fraction evolution

Compare predicted $f_{\text{tw}}$ against:

- EBSD twin boundary analysis,
- TEM observations,
- optical microscopy.

### 31.6 Strain-rate sensitivity

Run simulations at different strain rates and verify that the power-law slip rule produces realistic rate dependence.

---

## 32. Possible Extensions

The model can be extended in several directions.

### 32.1 Latent hardening

Introduce a hardening matrix:

$$
\dot{g}^{\alpha}
=
\sum_{\beta}
h^{\alpha\beta}
|\dot{\gamma}^{\beta}|
$$

This improves multislip and texture predictions.

### 32.2 Cubic elastic anisotropy

Replace isotropic Hooke elasticity with cubic elasticity using:

$$
C_{11}, C_{12}, C_{44}
$$

### 32.3 Full finite-strain hyperelastic formulation

Implement a multiplicative decomposition with an explicit update of:

$$
\mathbf{F}^{e}, \mathbf{F}^{p}, \mathbf{F}^{\text{tr}}, \mathbf{F}^{\text{tw}}
$$

### 32.4 Stress-state-dependent transformation kinetics

Introduce triaxiality and Lode angle dependence:

$$
K_J = K_J(\eta, \bar{\theta})
$$

or

$$
\alpha_{\text{OC}} =
\alpha_{\text{OC}}(\eta, \bar{\theta})
$$

### 32.5 Thermomechanical coupling

Update temperature from plastic work:

$$
\Delta T
=
\frac{\beta_{\text{Taylor}}}
{\rho c_p}
\boldsymbol{\sigma} :
\Delta \boldsymbol{\varepsilon}^{p}
$$

where $\beta_{\text{Taylor}}$ is the Taylor–Quinney coefficient.

### 32.6 Two-phase crystal plasticity

Represent martensite as a separate phase with its own slip systems:

- FCC austenite slip,
- BCT/BCC martensite slip,
- phase-specific hardening,
- load partitioning.

### 32.7 Transformation variant selection

Introduce martensite variants and transformation shear tensors:

$$
\boldsymbol{\varepsilon}^{\text{tr}}
=
\sum_{v}
f^{v}
\boldsymbol{\varepsilon}^{\text{tr},v}
$$

### 32.8 Twin reorientation

Implement predominant twin reorientation, PTR, or crystallographic reorientation of twinned regions.

### 32.9 Damage and fracture

Couple the model with continuum damage mechanics to predict ductile fracture after extensive transformation and twinning.

---

## 33. Model Selection Summary

| Model component | Selected model | Main reason |
|---|---|---|
| Elasticity | Isotropic Hooke | Efficient, stable, sufficient for small elastic strains |
| Plastic flow | Crystal plasticity slip systems | Captures crystallographic mechanisms directly |
| Slip systems | FCC `{111}<110>` | Standard austenite slip family |
| Slip kinetics | Power-law viscoplasticity | Rate dependence, explicit stability |
| Slip hardening | Voce per system | Captures saturation hardening |
| Macroscopic yield surface | Not used | Yielding emerges from slip-system CRSS |
| Martensite kinetics | Olson–Cohen + JMAK | Combines shear-band nucleation and sigmoidal growth |
| Transformation strain | Isotropic volumetric eigenstrain | Captures dominant volume expansion robustly |
| Twinning systems | FCC `{111}<112>` | Standard FCC twin family |
| Twinning criterion | Polar CRSS | Captures directional nature of twinning |
| Twinning kinetics | TANH law | Smooth, stable, saturating |
| Twin–martensite coupling | Linear coupling coefficient | Captures twin-assisted martensite nucleation |
| Numerical stabilization | Caps and limiters | Required for robust Abaqus/Explicit simulations |

---

## 34. Practical Usage Recommendations

1. **Use consistent units.**  
   If stress is in MPa, density and geometry must be consistent with MPa.

2. **Start with a single-element test.**  
   Verify stability before using the model in a large RVE.

3. **Use small stable time increments.**  
   Explicit simulations are sensitive to rapid internal variable changes.

4. **Monitor phase fractions.**  
   Ensure $f_M$ and $f_{\text{tw}}$ remain physically reasonable.

5. **Check orientation assignment.**  
   Incorrect orientations can produce unrealistic texture and yield behavior.

6. **Calibrate transformation and twinning separately if possible.**  
   Use phase fraction data for transformation and twin fraction data for twinning.

7. **Disable coupling initially.**  
   Set `alpha_TM = 0` first, then enable it after the base model is calibrated.

---

## 35. Known Items to Verify Before Publication

Before using this subroutine for publication or industrial prediction, verify:

1. Initial yield stress against Schmid-factor prediction.
2. Hardening curve against experimental tensile data.
3. Martensite fraction evolution against XRD/EBSD data.
4. Twin fraction evolution against microstructure measurements.
5. Strain-rate sensitivity against split-Hopkinson or multi-rate tests.
6. Sensitivity to explicit time increment and mass scaling.
7. Orientation dependence in single-crystal benchmark tests.
8. Volume expansion effect under constrained boundary conditions.
9. Energy outputs if energy balance is important.
10. Twin-system sign convention relative to the intended crystallographic texture.

---

## 36. Summary

This VUMAT provides a robust, mechanism-based constitutive model for TRIP/TWIP steels in Abaqus/Explicit. It combines:

- rate-dependent FCC crystal plasticity,
- Voce hardening,
- Olson–Cohen shear-band transformation nucleation,
- JMAK strain-driven transformation kinetics,
- isotropic transformation volumetric expansion,
- polar CRSS twinning activation,
- smooth TANH twin kinetics,
- twin-assisted martensite nucleation,
- explicit numerical stabilization.

The model is not intended to be a fully resolved two-phase micromechanical model. Instead, it is a computationally efficient, physically motivated material-point model capable of capturing the essential coupling between slip, twinning, and martensitic transformation in metastable austenitic/TRIP steels.

It is especially suitable for simulations where the simultaneous evolution of plasticity, phase transformation, and twinning controls the macroscopic mechanical response.

---

## 37. Key References for Further Reading

1. J. R. Rice, crystal plasticity foundations.
2. D. Peirce, R. J. Asaro, A. Needleman, rate-dependent crystal plasticity.
3. S. R. Kalidindi, C. A. Bronkhorst, L. Anand, crystal plasticity calibration for FCC metals.
4. Y. Huang, Harvard crystal plasticity VUMAT/UMAT implementation.
5. G. B. Olson, M. Cohen, strain-induced martensitic transformation.
6. D. P. Koistinen, R. E. Marburger, martensite transformation kinetics.
7. W. A. Johnson, R. F. Mehl, M. Avrami, A. N. Kolmogorov, JMAK transformation theory.
8. Y. Tomita, T. Iwamoto, TRIP steel constitutive modeling.
9. P. Van Houtte, predominant twin reorientation.
10. Advanced high-strength steel literature for TRIP/TWIP mechanisms.

---

## 38. License and Usage Note

This model is intended for research and educational use. Industrial application requires careful calibration and validation against experimental data for the specific steel grade, heat treatment, texture, strain rate, and temperature range.