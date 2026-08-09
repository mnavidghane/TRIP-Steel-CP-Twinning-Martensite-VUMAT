Coupled Crystal Plasticity–Martensitic Transformation–Twinning VUMAT for TRIP Steels
Abaqus/Explicit user material subroutine for transformation-induced plasticity, mechanical twinning, and crystallographic slip in metastable austenitic/TRIP steels

1. Executive Summary
This repository contains a VUMAT user material subroutine for Abaqus/Explicit that implements a coupled constitutive model for TRIP steels and related metastable austenitic alloys. The model is designed to capture the simultaneous contribution of:
Crystallographic slip in FCC austenite using a rate-dependent crystal plasticity formulation.
Martensitic transformation from metastable retained austenite to martensite during plastic deformation.
Mechanical twinning on FCC twin systems.
Twin–martensite interaction, where twin formation can accelerate martensite nucleation.
Transformation-induced volumetric expansion, represented through an isotropic transformation eigenstrain.
Numerical stabilization features required for robust explicit finite element simulations.
The subroutine is written for a homogenized material-point description. Austenite is treated as the crystallographically active FCC phase carrying slip and twinning. Martensite is represented primarily by its volume fraction, its volumetric transformation strain, and its indirect contribution to the macroscopic hardening response. A full two-phase FCC/BCT crystal plasticity treatment is not implemented in the current version but is a natural extension.
The model is suitable for:
TRIP steel forming simulations.
Crashworthiness and energy-absorption studies.
Polycrystalline RVE simulations when combined with grain orientations.
Investigations of phase fraction evolution under mechanical loading.
ICME workflows linking microstructure evolution to macroscopic mechanical response.
Important modeling note:
The theoretical framework behind this work can be formulated in a finite-strain multiplicative kinematic setting. The present VUMAT uses an incrementally additive, corotational/hypoelastic-style implementation compatible with Abaqus/Explicit. It is robust and efficient for engineering simulations where elastic strains remain small compared with plastic, transformation, and twinning strains.
2. Physical Mechanisms Represented
Physical mechanism
Representation in the VUMAT
Primary state variables
Elastic response
Isotropic Hookean elasticity
Stress tensor
Crystallographic slip
12 FCC {111}<110> slip systems, rate-dependent power-law flow, Voce hardening
Slip shear strains, CRSS values
Martensitic transformation
Olson–Cohen-type shear-band nucleation + JMAK-type strain-driven kinetics
Martensite volume fraction f_M
Transformation eigenstrain
Isotropic volumetric expansion proportional to Δf_M
Accumulated transformation strain eps_M
Mechanical twinning
12 FCC {111}<112> twin systems, polar CRSS activation, TANH kinetics
Twin volume fraction f_tw
Twin-induced plastic strain
Twin shear distributed over active twin systems
Twin plastic strain tensor
Twin–martensite interaction
Additional martensite nucleation proportional to twin fraction increment
Coupling parameter alpha_TM
Phase fraction constraint
f_M + f_tw <= 1 enforced
f_M, f_tw
Numerical stabilization
Slip increment cap, volume-fraction rate caps, anti-overshoot limiter
Internal algorithmic variables
3. Constitutive Architecture
The model belongs to the class of mechanism-based crystal plasticity models with internal phase evolution. Instead of using a macroscopic yield surface such as von Mises, Hill48, or Barlat, the onset of plastic flow is governed locally on each crystallographic system by the resolved shear stress and the corresponding critical resolved shear stress.
The overall inelastic strain rate is conceptually decomposed as
ε
˙
inel
=
ε
˙
slip
+
ε
˙
twin
+
ε
˙
tr
ε
˙
  
inel
 = 
ε
˙
  
slip
 + 
ε
˙
  
twin
 + 
ε
˙
  
tr
 
where
ε
˙
slip
ε
˙
  
slip
  is the plastic strain rate from dislocation glide,
ε
˙
twin
ε
˙
  
twin
  is the pseudo-plastic strain rate associated with mechanical twinning,
ε
˙
tr
ε
˙
  
tr
  is the transformation strain rate associated with martensite formation.
In the current implementation, the strain increment passed by Abaqus is treated incrementally:
Δ
ε
=
Δ
ε
e
+
Δ
ε
slip
+
Δ
ε
twin
+
Δ
ε
tr,vol
Δε=Δε 
e
 +Δε 
slip
 +Δε 
twin
 +Δε 
tr,vol
 
The stress update is then performed from the elastic part only:
Δ
σ
=
D
e
:
(
Δ
ε
−
Δ
ε
slip
−
Δ
ε
twin
−
Δ
ε
tr,vol
)
Δσ=D 
e
 :(Δε−Δε 
slip
 −Δε 
twin
 −Δε 
tr,vol
 )
where 
D
e
D 
e
  is the isotropic elastic stiffness tensor.
4. Kinematic Framework
4.1 Conceptual finite-strain picture
A rigorous finite-strain description of TRIP/TWIP behavior often uses a multiplicative decomposition of the deformation gradient:
F
=
F
e
F
p
F
tr
F
tw
F=F 
e
 F 
p
 F 
tr
 F 
tw
 
where
F
e
F 
e
  is the elastic lattice distortion,
F
p
F 
p
  is the plastic deformation due to crystallographic slip,
F
tr
F 
tr
  is the transformation deformation due to martensite formation,
F
tw
F 
tw
  is the twinning deformation.
This decomposition is useful because it separates physically distinct mechanisms:
Component
Mechanism
F
e
F 
e
 
Elastic stretching and lattice rotation
F
p
F 
p
 
Dislocation slip
F
tr
F 
tr
 
Martensitic transformation strain
F
tw
F 
tw
 
Twinning shear and reorientation
4.2 Implemented incremental form
The present VUMAT does not explicitly track all deformation gradients. Instead, it uses an incremental strain decomposition in the corotational frame supplied by Abaqus/Explicit:
Δ
ε
=
Δ
ε
e
+
Δ
ε
cp
+
Δ
ε
tw
+
Δ
ε
tr,vol
Δε=Δε 
e
 +Δε 
cp
 +Δε 
tw
 +Δε 
tr,vol
 
This approach is computationally efficient and robust for explicit simulations. It is particularly appropriate when:
elastic strains are small,
plastic strains are large,
transformation strains are moderate,
the main interest is macroscopic stress–strain response and phase fraction evolution.
For problems requiring exact finite-strain hyperelastic consistency, strong elastic anisotropy, or explicit lattice reorientation, a full multiplicative implementation would be preferable.
5. Elastic Response
5.1 Implemented model
The subroutine uses isotropic linear elasticity:
λ
=
E
ν
(
1
+
ν
)
(
1
−
2
ν
)
λ= 
(1+ν)(1−2ν)
Eν
​
 
μ
=
E
2
(
1
+
ν
)
μ= 
2(1+ν)
E
​
 
K
=
λ
+
2
3
μ
K=λ+ 
3
2
​
 μ
The elastic stiffness matrix in Voigt notation is assembled as:
D
e
=
[
λ
+
2
μ
λ
λ
0
0
0
λ
λ
+
2
μ
λ
0
0
0
λ
λ
λ
+
2
μ
0
0
0
0
0
0
2
μ
0
0
0
0
0
0
2
μ
0
0
0
0
0
0
2
μ
]
D 
e
 = 
​
  
λ+2μ
λ
λ
0
0
0
​
  
λ
λ+2μ
λ
0
0
0
​
  
λ
λ
λ+2μ
0
0
0
​
  
0
0
0
2μ
0
0
​
  
0
0
0
0
2μ
0
​
  
0
0
0
0
0
2μ
​
  
​
 
The stress increment is
Δ
σ
=
D
e
:
Δ
ε
e
Δσ=D 
e
 :Δε 
e
 
5.2 Elasticity model alternatives
Elastic model
Advantages
Disadvantages
Suitability for present code
Isotropic Hooke
Simple, robust, only two parameters
Ignores crystal elastic anisotropy
Selected
Saint Venant–Kirchhoff
Finite-strain compatible
Poor for large strains, not ideal for metals
Not selected
Neo-Hookean hyperelastic
Thermodynamically consistent for finite elastic strains
Still isotropic, more complex
Possible extension
Cubic anisotropic elasticity
Captures FCC elastic anisotropy
Requires 
C
11
,
C
12
,
C
44
C 
11
​
 ,C 
12
​
 ,C 
44
​
 
Future extension
Anisotropic hyperelasticity
General texture-dependent elasticity
Many parameters, calibration burden
Not selected
5.3 Why isotropic elasticity was selected
The current model focuses on the dominant sources of inelasticity in TRIP steels:
dislocation slip,
martensitic transformation,
mechanical twinning,
transformation/twinning-induced hardening.
For many engineering forming and impact simulations, the elastic strain is small compared with the plastic and transformation strains. Therefore, isotropic elasticity provides a good compromise between:
numerical stability,
computational efficiency,
ease of calibration,
robustness in explicit simulations.
The main plastic anisotropy is captured through the crystallographic slip and twin systems, not through elastic anisotropy.
6. Plastic Flow: Why Crystal Plasticity Instead of a Macroscopic Yield Criterion?
Several macroscopic yield criteria could be used to describe plastic yielding in steels:
Yield/plasticity model
Strengths
Weaknesses
von Mises 
J
2
J 
2
​
 
Simple, robust, isotropic
Cannot capture crystallographic texture or orientation effects
Tresca
Simple shear-based criterion
Corners in yield surface, less accurate for FCC metals
Drucker–Prager
Includes hydrostatic pressure sensitivity
Mainly for soils, concrete, porous materials
Hill48
Captures orthotropic sheet anisotropy
Empirical, no direct slip-system information
Barlat
Advanced sheet anisotropy
More parameters, higher complexity
Crystal plasticity
Directly models slip systems, texture, orientation, CRSS evolution
Higher computational cost
6.1 Selected model: rate-dependent crystal plasticity
For TRIP steels, the physically dominant plastic mechanism in austenite is slip on FCC systems:
{
111
}
⟨
110
⟩
{111}⟨110⟩
There are 12 such slip systems. Crystal plasticity is selected because it can represent:
crystallographic orientation,
texture evolution tendencies,
grain-level anisotropy,
slip-system-level hardening,
coupling with twinning,
coupling with phase transformation,
localized deformation in polycrystalline RVEs.
In a crystal plasticity framework, the macroscopic yield surface is not imposed directly. Instead, yielding emerges from the collective activation of slip systems according to Schmid-type resolved shear stress conditions.
7. Crystallographic Slip Model
7.1 Slip systems
The model uses the standard FCC slip family:
{
111
}
⟨
110
⟩
{111}⟨110⟩
with 12 systems:
Plane family
Slip directions
Number of systems
(
111
)
(111)
[
0
1
ˉ
1
]
[0 
1
ˉ
 1], 
[
10
1
ˉ
]
[10 
1
ˉ
 ], 
[
1
ˉ
10
]
[ 
1
ˉ
 10]
3
(
1
1
ˉ
1
)
(1 
1
ˉ
 1)
[
011
]
[011], 
[
101
]
[101], 
[
1
ˉ
1
ˉ
0
]
[ 
1
ˉ
  
1
ˉ
 0]
3
(
1
ˉ
11
)
( 
1
ˉ
 11)
[
01
1
ˉ
]
[01 
1
ˉ
 ], 
[
101
]
[101], 
[
1
ˉ
10
]
[ 
1
ˉ
 10]
3
(
11
1
ˉ
)
(11 
1
ˉ
 )
[
011
]
[011], 
[
10
1
ˉ
]
[10 
1
ˉ
 ], 
[
1
ˉ
1
ˉ
0
]
[ 
1
ˉ
  
1
ˉ
 0]
3
The slip systems are hard-coded in the local material coordinate system. For polycrystalline simulations, grain orientation must be supplied through Abaqus *ORIENTATION.
7.2 Schmid tensor
For each slip system 
α
α, with slip direction 
s
α
s 
α
  and slip-plane normal 
n
α
n 
α
 , the symmetric Schmid tensor is
M
α
=
1
2
(
s
α
⊗
n
α
+
n
α
⊗
s
α
)
M 
α
 = 
2
1
​
 (s 
α
 ⊗n 
α
 +n 
α
 ⊗s 
α
 )
The resolved shear stress is
τ
α
=
M
α
:
σ
τ 
α
 =M 
α
 :σ
7.3 Rate-dependent power-law flow rule
The shear rate on slip system 
α
α is
γ
˙
α
=
γ
˙
0
∣
τ
α
g
α
∣
1
/
m
sign
⁡
(
τ
α
)
γ
˙
​
  
α
 = 
γ
˙
​
  
0
​
  
​
  
g 
α
 
τ 
α
 
​
  
​
  
1/m
 sign(τ 
α
 )
where
γ
˙
0
γ
˙
​
  
0
​
  is the reference shear strain rate,
m
m is the rate-sensitivity exponent,
g
α
g 
α
  is the current critical resolved shear stress, CRSS.
The slip increment is
Δ
γ
α
=
γ
˙
α
Δ
t
Δγ 
α
 = 
γ
˙
​
  
α
 Δt
7.4 Plastic strain from slip
The plastic strain increment due to slip is
Δ
ε
cp
=
∑
α
=
1
12
M
α
Δ
γ
α
Δε 
cp
 = 
α=1
∑
12
​
 M 
α
 Δγ 
α
 
7.5 Voce hardening law
The CRSS of each slip system evolves according to a Voce-type saturation law:
g
α
(
Γ
α
)
=
τ
0
+
(
τ
s
−
τ
0
)
[
1
−
exp
⁡
(
−
h
0
Γ
α
τ
s
−
τ
0
)
]
g 
α
 (Γ 
α
 )=τ 
0
​
 +(τ 
s
​
 −τ 
0
​
 )[1−exp(− 
τ 
s
​
 −τ 
0
​
 
h 
0
​
 Γ 
α
 
​
 )]
where
τ
0
τ 
0
​
  is the initial CRSS,
τ
s
τ 
s
​
  is the saturation CRSS,
h
0
h 
0
​
  is the initial hardening rate,
Γ
α
Γ 
α
  is the accumulated absolute slip on system 
α
α.
If 
τ
s
≤
τ
0
τ 
s
​
 ≤τ 
0
​
 , the model falls back to linear hardening:
g
α
=
τ
0
+
h
0
Γ
α
g 
α
 =τ 
0
​
 +h 
0
​
 Γ 
α
 
7.6 Hardening model alternatives
Hardening model
Advantages
Disadvantages
Reason for selection/rejection
Linear hardening
Simple
No saturation
Fallback only
Power-law hardening
Common for metals
Less natural saturation behavior
Not selected
Voce hardening
Captures saturation typical of FCC metals
Requires 
τ
s
,
h
0
τ 
s
​
 ,h 
0
​
 
Selected
Isotropic macroscopic hardening
Simple continuum model
No slip-system resolution
Not suitable for CP
Kinematic hardening
Captures Bauschinger effect
More parameters, not system-based here
Not implemented
Latent hardening matrix
Captures cross-system interaction
Requires latent hardening coefficients
Future extension
7.7 Important slip-hardening simplification
The current implementation uses self-hardening only. The CRSS of each slip system evolves according to its own accumulated slip:
Γ
α
=
∑
∣
Δ
γ
α
∣
Γ 
α
 =∑∣Δγ 
α
 ∣
A more advanced model would include a latent hardening matrix:
g
˙
α
=
∑
β
h
α
β
∣
γ
˙
β
∣
g
˙
​
  
α
 = 
β
∑
​
 h 
αβ
 ∣ 
γ
˙
​
  
β
 ∣
Latent hardening is often important for texture evolution and multislip behavior in FCC metals. It is a recommended extension.
8. Martensitic Transformation
Martensitic transformation is one of the central mechanisms in TRIP steels. Metastable retained austenite transforms into martensite during plastic deformation, producing:
additional hardening,
delayed necking,
increased strength,
improved energy absorption,
local volumetric expansion,
transformation-induced plastic strain.
The transformation model in this VUMAT is built from three parts:
Transformation kinetics,
Transformation eigenstrain,
Coupling with twinning and plastic strain.
9. Thermodynamic Background of Martensitic Transformation
Although the present VUMAT uses a robust strain-driven implementation, the physical basis of martensitic transformation is thermodynamic.
9.1 Chemical driving force
The chemical free-energy difference between austenite and martensite can be written approximately as
Δ
G
chem
=
Δ
H
M
−
T
Δ
S
M
ΔG 
chem
​
 =ΔH 
M
​
 −TΔS 
M
​
 
where
Δ
H
M
ΔH 
M
​
  is the enthalpy change,
Δ
S
M
ΔS 
M
​
  is the entropy change,
T
T is absolute temperature.
Transformation becomes favorable when the total driving force exceeds a critical energy barrier.
9.2 Mechanical driving force
An applied stress state can assist transformation by performing mechanical work on the transformation strain:
Δ
G
mech
=
σ
:
ε
tr
ΔG 
mech
​
 =σ:ε 
tr
 
Thus the total driving force may be written conceptually as
Δ
G
M
=
Δ
G
chem
+
Δ
G
mech
ΔG 
M
​
 =ΔG 
chem
​
 +ΔG 
mech
​
 
or, in a more explicit form,
Δ
G
M
=
Δ
H
M
−
T
Δ
S
M
+
σ
:
ε
tr
ΔG 
M
​
 =ΔH 
M
​
 −TΔS 
M
​
 +σ:ε 
tr
 
9.3 Critical barrier and Olson–Cohen-type condition
In Olson–Cohen-type physical nucleation models, transformation starts when the driving force exceeds a critical barrier:
Δ
G
M
≥
G
M
crit
ΔG 
M
​
 ≥G 
M
crit
​
 
A common simplified form is
G
M
crit
=
G
M
0
−
c
σ
2
G 
M
crit
​
 =G 
M
0
​
 −cσ 
2
 
where
G
M
0
G 
M
0
​
  is the stress-free barrier,
c
c is a material constant,
σ
σ represents the applied stress intensity.
This expression captures the idea that mechanical loading lowers the nucleation barrier.
9.4 Stress-assisted versus strain-induced transformation
Two regimes are commonly distinguished:
Regime
Description
Stress-assisted transformation
Transformation is assisted by applied stress before large plastic strain; often important near or below 
M
s
M 
s
​
 
Strain-induced transformation
Transformation is driven by plastic deformation, shear bands, and defect generation; dominant in many TRIP steels
The present VUMAT primarily captures the strain-induced transformation regime through accumulated plastic strain and shear-band-like nucleation variables.
9.5 Stress triaxiality and Lode angle
Advanced transformation models may include the effect of:
stress triaxiality,
Lode angle,
hydrostatic pressure,
shear-dominant versus tension-dominant loading.
Stress triaxiality is often defined as
η
=
σ
m
σ
ˉ
η= 
σ
ˉ
 
σ 
m
​
 
​
 
where
σ
m
=
1
3
tr
⁡
(
σ
)
σ 
m
​
 = 
3
1
​
 tr(σ)
and 
σ
ˉ
σ
ˉ
  is the von Mises equivalent stress.
High positive triaxiality generally promotes volumetric expansion associated with martensitic transformation, while shear-dominant states may alter the transformation rate and variant selection.
The current VUMAT does not explicitly compute triaxiality or Lode angle inside the kinetic law. However, these effects can be introduced in future versions by modifying the transformation kinetics parameters as functions of the local stress state.
9.6 Temperature and adiabatic heating
Martensitic transformation is temperature dependent. Increasing temperature generally stabilizes austenite and reduces the transformation rate. At high strain rates, adiabatic heating may raise the local temperature and suppress transformation.
The current implementation stores temperature as a state variable but keeps it fixed by default. Temperature dependence can be introduced externally by making parameters such as:
K_J,
tau_twin,
tau0,
taus,
beta_OC,
functions of temperature.
10. Martensitic Transformation Kinetics
10.1 General kinetic behavior
The evolution of martensite volume fraction versus plastic strain in TRIP steels is often sigmoidal:
an initial incubation stage,
a rapid transformation stage,
a saturation stage as austenite is consumed.
A robust kinetic model should capture these three stages.
The present model combines two kinetic descriptions:
Olson–Cohen-type shear-band nucleation,
JMAK-type transformation kinetics.
They are combined as independent transformation mechanisms.
11. Olson–Cohen-Type Shear-Band Nucleation Model
The Olson–Cohen model is one of the most physically meaningful models for strain-induced martensitic transformation. It assumes that plastic deformation generates shear bands, and intersections of shear bands act as martensite nucleation sites.
11.1 Shear-band density evolution
The VUMAT uses a normalized shear-band density variable 
N
sb
N 
sb
​
 . Its increment is
Δ
N
sb
=
α
OC
(
1
−
N
sb
)
Δ
ε
ˉ
p
ε
0
,
OC
ΔN 
sb
​
 =α 
OC
​
 (1−N 
sb
​
 ) 
ε 
0,OC
​
 
Δ 
ε
ˉ
  
p
 
​
 
where
α
OC
α 
OC
​
  controls the shear-band generation rate,
ε
0
,
OC
ε 
0,OC
​
  is a reference plastic strain,
Δ
ε
ˉ
p
Δ 
ε
ˉ
  
p
  is the equivalent plastic strain increment.
The updated shear-band density is capped between 0 and 1:
0
≤
N
sb
≤
1
0≤N 
sb
​
 ≤1
11.2 Martensite fraction from Olson–Cohen nucleation
The transformed fraction associated with shear-band nucleation is
f
M
OC
=
1
−
exp
⁡
[
−
β
OC
N
sb
n
OC
]
f 
M
OC
​
 =1−exp[−β 
OC
​
 N 
sb
n 
OC
​
 
​
 ]
where
β
OC
β 
OC
​
  controls nucleation efficiency,
n
OC
n 
OC
​
  controls the shape of the transformation curve.
11.3 Advantages of the Olson–Cohen model
Physically based on shear-band nucleation.
Captures strain-induced transformation.
Naturally produces sigmoidal behavior.
Compatible with crystal plasticity.
Parameters have microstructural interpretation.
11.4 Limitations
Does not explicitly include stress triaxiality or Lode angle.
Does not explicitly include temperature unless parameters are modified.
Does not distinguish individual martensite variants.
Assumes austenite transformation is controlled primarily by plastic strain and shear-band density.
12. JMAK-Type Transformation Kinetics
The Johnson–Mehl–Avrami–Kolmogorov model is a classical nucleation-and-growth model. In strain-driven form, it may be written as
f
M
JMAK
=
1
−
exp
⁡
[
−
K
J
(
ε
ˉ
p
)
n
J
]
f 
M
JMAK
​
 =1−exp[−K 
J
​
 ( 
ε
ˉ
  
p
 ) 
n 
J
​
 
 ]
where
K
J
K 
J
​
  is a temperature-dependent rate constant,
n
J
n 
J
​
  is the Avrami exponent,
ε
ˉ
p
ε
ˉ
  
p
  is the accumulated equivalent plastic strain.
12.1 Advantages of JMAK
Captures sigmoidal transformation behavior.
Simple and robust.
Useful for fitting experimental transformation curves.
Can absorb temperature dependence through 
K
J
K 
J
​
 .
12.2 Limitations
Less mechanistic than Olson–Cohen.
Does not explicitly represent shear-band nucleation.
Does not directly include stress-state effects.
Requires calibration against phase fraction data.
13. Combined Olson–Cohen and JMAK Kinetics
The VUMAT combines the two transformation descriptions as independent mechanisms. If two independent transformation mechanisms produce fractions 
f
M
OC
f 
M
OC
​
  and 
f
M
JMAK
f 
M
JMAK
​
 , the combined fraction is
f
M
base
=
1
−
(
1
−
f
M
OC
)
(
1
−
f
M
JMAK
)
f 
M
base
​
 =1−(1−f 
M
OC
​
 )(1−f 
M
JMAK
​
 )
This expression means that the untransformed austenite fraction is the product of the untransformed fractions from each mechanism:
1
−
f
M
base
=
(
1
−
f
M
OC
)
(
1
−
f
M
JMAK
)
1−f 
M
base
​
 =(1−f 
M
OC
​
 )(1−f 
M
JMAK
​
 )
13.1 Why combine OC and JMAK?
Model
Captures
Olson–Cohen
Shear-band nucleation, strain-induced transformation
JMAK
Overall sigmoidal nucleation/growth behavior
Combined model
More flexible representation of experimental data
The combined model allows the user to represent materials where transformation is controlled by both defect generation and overall nucleation/growth statistics.
14. Twin–Martensite Coupling
Mechanical twins can act as additional nucleation sites for martensite. To represent this, the model adds a coupling term:
f
M
base
←
f
M
base
+
α
T
M
Δ
f
tw
(
1
−
f
M
base
)
f 
M
base
​
 ←f 
M
base
​
 +α 
TM
​
 Δf 
tw
​
 (1−f 
M
base
​
 )
where
α
T
M
α 
TM
​
  is the twin-to-martensite coupling coefficient,
Δ
f
tw
Δf 
tw
​
  is the twin volume fraction increment,
1
−
f
M
base
1−f 
M
base
​
  is the remaining transformable austenite.
If 
α
T
M
=
0
α 
TM
​
 =0, the coupling is disabled.
14.1 Physical meaning
This term captures the idea that twin boundaries can:
act as martensite nucleation sites,
increase local defect density,
assist transformation in metastable austenite,
alter local transformation kinetics.
14.2 Monotonicity constraint
Martensite transformation is treated as irreversible:
f
M
new
≥
f
M
old
f 
M
new
​
 ≥f 
M
old
​
 
This prevents nonphysical reversal of martensite to austenite during unloading.
15. Transformation Eigenstrain
Martensitic transformation in steels is accompanied by a volume expansion. The present model represents this using an isotropic volumetric eigenstrain.
Let
Δ
V
/
V
ΔV/V
be the volumetric transformation strain associated with complete transformation. In the code this parameter is dVoV. Typical values for austenite-to-martensite transformation are approximately
0.02
≤
Δ
V
/
V
≤
0.04
0.02≤ΔV/V≤0.04
For a martensite fraction increment 
Δ
f
M
Δf 
M
​
 , the volumetric eigenstrain increment per spatial direction is
Δ
ε
i
i
tr
=
Δ
f
M
Δ
V
/
V
3
Δε 
ii
tr
​
 = 
3
Δf 
M
​
 ΔV/V
​
 
The corresponding stress correction for constrained expansion is
Δ
σ
i
i
=
−
K
Δ
f
M
Δ
V
/
V
Δσ 
ii
​
 =−KΔf 
M
​
 ΔV/V
for the three normal components.
15.1 Transformation effect model alternatives
Transformation strain model
Description
Advantages
Disadvantages
Unity model
No transformation strain
Simplest
Physically inadequate for TRIP steels
Approximate isotropic volumetric model
Small-strain volumetric expansion
Simple
Less exact for finite strain
Exact volumetric model
Uses transformation Jacobian 
J
tr
J 
tr
 
Thermodynamically cleaner
Still isotropic
Directional shear model
Includes variant-specific shear
Captures variant selection
Requires variant data and calibration
Combined volume-shear model
Includes dilatation and shear
Most complete
Complex, many parameters
Statistical/microstructural model
Spatially heterogeneous transformation
Captures local effects
High computational cost
15.2 Why the volumetric eigenstrain model was selected
The current model uses an isotropic volumetric transformation strain because:
The dominant macroscopic effect of austenite-to-martensite transformation is volume expansion.
The model remains robust for explicit simulations.
It requires only one main parameter, dVoV.
It avoids the need to track individual martensite variants.
It can be calibrated from dilatometry or crystallographic data.
A directional shear or combined volume-shear transformation model is a recommended extension when variant selection, transformation texture, or shear-assisted transformation is important.
16. Mechanical Twinning
Mechanical twinning is another important deformation mechanism in low-SFE FCC alloys, TWIP steels, and some TRIP-assisted steels. Twinning contributes to:
plastic strain,
work hardening,
texture evolution,
dynamic Hall–Petch-type strengthening,
martensite nucleation at twin boundaries.
The present model represents twinning using 12 FCC twin systems:
{
111
}
⟨
112
⟩
{111}⟨112⟩
17. Twinning Kinematic Model
17.1 Twin Schmid tensors
For each twin system 
t
t, a symmetric Schmid-like tensor is constructed:
M
t
=
1
2
(
s
t
⊗
n
t
+
n
t
⊗
s
t
)
M 
t
 = 
2
1
​
 (s 
t
 ⊗n 
t
 +n 
t
 ⊗s 
t
 )
where
n
t
n 
t
  is the twin-plane normal,
s
t
s 
t
  is the twin shear direction.
The resolved twin shear stress is
τ
t
=
M
t
:
σ
τ 
t
 =M 
t
 :σ
17.2 Polar CRSS activation
Twinning is inherently directional. Unlike ordinary slip, which can often be treated as approximately bidirectional, twinning occurs in a specific crystallographic sense.
The model therefore uses a polar CRSS criterion:
τ
t
>
τ
c
tw
τ 
t
 >τ 
c
tw
​
 
and also requires
τ
t
>
0
τ 
t
 >0
Thus a twin system activates only when the resolved shear stress is both positive and larger than the critical twin shear stress.
This is more physical than a non-polar absolute-value criterion.
17.3 Twinning condition model alternatives
Twinning criterion
Advantages
Disadvantages
Absolute CRSS
Simple
Does not distinguish twin sense
Polar CRSS
Captures directional twinning
Requires sign convention care
Energy-based SFE criterion
Thermodynamically richer
Requires stacking fault energy data
Full nucleation theory
Most physical
Complex and data-intensive
The selected model is the polar CRSS criterion because it is simple, robust, and physically adequate for engineering-scale simulations.
18. Twinning Kinetics
18.1 TANH kinetic law
The model uses a smooth hyperbolic tangent law to determine the target twin volume fraction:
f
tw
target
=
A
tw
1
2
[
1
+
tanh
⁡
(
τ
max
⁡
t
−
τ
c
tw
B
tw
)
]
f 
tw
target
​
 =A 
tw
​
  
2
1
​
 [1+tanh( 
B 
tw
​
 
τ 
max
t
​
 −τ 
c
tw
​
 
​
 )]
where
A
tw
A 
tw
​
  is the maximum twin volume fraction,
B
tw
B 
tw
​
  controls the transition width,
τ
max
⁡
t
τ 
max
t
​
  is the maximum resolved twin shear stress among all twin systems,
τ
c
tw
τ 
c
tw
​
  is the critical twin shear stress.
18.2 Irreversibility
Twin growth is treated as irreversible:
Δ
f
tw
=
max
⁡
(
0
,
f
tw
target
−
f
tw
old
)
Δf 
tw
​
 =max(0,f 
tw
target
​
 −f 
tw
old
​
 )
This prevents detwinning during unloading in the current implementation.
18.3 Rate cap for numerical stability
The twin fraction increment is limited by a maximum value per increment:
Δ
f
tw
≤
Δ
f
tw
max
⁡
Δf 
tw
​
 ≤Δf 
tw
max
​
 
In the current code, the default cap is
Δ
f
tw
max
⁡
=
0.02
Δf 
tw
max
​
 =0.02
This improves explicit stability and prevents abrupt strain spikes.
18.4 Twin plastic strain
The twin plastic strain increment is distributed over active twin systems according to their relative driving stresses.
Define the active-system weight:
w
t
=
max
⁡
(
0
,
τ
t
−
τ
c
tw
)
w 
t
 =max(0,τ 
t
 −τ 
c
tw
​
 )
The total weight is
W
=
∑
t
w
t
W= 
t
∑
​
 w 
t
 
If 
W
>
0
W>0, the twin plastic strain increment is
Δ
ε
tw
=
∑
t
w
t
W
Δ
f
tw
eff
γ
tw
M
t
Δε 
tw
 = 
t
∑
​
  
W
w 
t
 
​
 Δf 
tw
eff
​
 γ 
tw
​
 M 
t
 
where
γ
tw
γ 
tw
​
  is the characteristic twin shear,
for FCC twinning, 
γ
tw
≈
0.707
γ 
tw
​
 ≈0.707.
18.5 Twinning kinetics alternatives
Kinetic model
Advantages
Disadvantages
Reason for selection/rejection
Linear kinetic law
Simple
No saturation, unrealistic at large strain
Not selected
Logistic kinetic law
Captures saturation
Requires rate parameters
Possible alternative
PTR, predominant twin reorientation
Captures texture reorientation
Numerically discontinuous, complex
Future extension
TANH kinetic law
Smooth, stable, continuous, captures gradual activation and saturation
Phenomenological
Selected
Fully coupled twin-slip crystal plasticity
Most physical
High cost and complexity
Future extension
18.6 Why TANH kinetics were selected
The hyperbolic tangent law was selected because it provides:
smooth activation,
continuous derivatives,
stable explicit integration,
natural saturation,
easy calibration,
no abrupt jumps in volume fraction.
This is particularly important in Abaqus/Explicit, where sudden internal strain increments can destabilize the solution.
19. Phase Fraction Constraint
Because both martensite and twins consume austenite, the model enforces
f
M
+
f
tw
≤
1
f 
M
​
 +f 
tw
​
 ≤1
The implementation proceeds as follows:
A preliminary twin increment is computed and capped by the available austenite.
Martensite transformation is computed, including twin-assisted coupling.
The twin fraction is re-checked against the updated martensite fraction.
If necessary, the twin fraction is reduced so that the total transformed fraction does not exceed unity.
This prevents nonphysical phase fractions.
20. Complete Incremental Algorithm
For each material point and each increment, the VUMAT performs the following steps.
Step 0: Initialization
If stepTime = 0, Abaqus is performing an initial check. The subroutine performs a purely elastic update and does not activate plasticity, twinning, or transformation.
If totalTime = 0, state variables are initialized:
slip strains = 0,
CRSS = tau0,
martensite fraction = 0,
twin fraction = 0,
equivalent plastic strain = 0,
shear-band density = 0,
transformation strain = 0,
temperature = 298 K by default.
Step 1: Elastic predictor
A trial stress is computed:
σ
trial
=
σ
old
+
D
e
:
Δ
ε
σ 
trial
 =σ 
old
 +D 
e
 :Δε
Step 2: Slip-system resolved shear stresses
For each of the 12 slip systems:
τ
α
=
M
α
:
σ
trial
τ 
α
 =M 
α
 :σ 
trial
 
Step 3: Rate-dependent slip increments
Δ
γ
α
=
γ
˙
0
∣
τ
α
g
α
∣
1
/
m
sign
⁡
(
τ
α
)
Δ
t
Δγ 
α
 = 
γ
˙
​
  
0
​
  
​
  
g 
α
 
τ 
α
 
​
  
​
  
1/m
 sign(τ 
α
 )Δt
The slip increments are capped for stability.
Step 4: Assemble slip plastic strain
Δ
ε
cp
=
∑
α
=
1
12
M
α
Δ
γ
α
Δε 
cp
 = 
α=1
∑
12
​
 M 
α
 Δγ 
α
 
Step 5: Anti-overshoot limiter
The deviatoric part of the slip strain increment is compared with the deviatoric part of the imposed strain increment. If the slip strain increment is larger, it is scaled back:
Δ
ε
cp
←
s
Δ
ε
cp
Δε 
cp
 ←sΔε 
cp
 
with
s
=
∥
dev
⁡
Δ
ε
∥
∥
dev
⁡
Δ
ε
cp
∥
s= 
∥devΔε 
cp
 ∥
∥devΔε∥
​
 
This prevents explicit instability caused by excessive plastic correction in a single increment.
Step 6: Update Voce hardening
For each slip system:
Γ
new
α
=
Γ
old
α
+
∣
Δ
γ
α
∣
Γ 
new
α
​
 =Γ 
old
α
​
 +∣Δγ 
α
 ∣
Then 
g
α
g 
α
  is updated using the Voce law.
Step 7: Twinning resolved shear stresses
For each twin system:
τ
t
=
M
t
:
σ
trial
τ 
t
 =M 
t
 :σ 
trial
 
Step 8: TANH twin target fraction
f
tw
target
=
A
tw
1
2
[
1
+
tanh
⁡
(
τ
max
⁡
t
−
τ
c
tw
B
tw
)
]
f 
tw
target
​
 =A 
tw
​
  
2
1
​
 [1+tanh( 
B 
tw
​
 
τ 
max
t
​
 −τ 
c
tw
​
 
​
 )]
Step 9: Twin fraction update
The twin fraction increment is made irreversible and capped.
Step 10: Martensite transformation
The equivalent plastic strain increment is computed from the slip strain increment:
Δ
ε
ˉ
p
=
2
3
Δ
ε
cp
:
Δ
ε
cp
Δ 
ε
ˉ
  
p
 = 
3
2
​
 Δε 
cp
 :Δε 
cp
 
​
 
Then the Olson–Cohen shear-band density is updated:
N
sb,new
=
N
sb,old
+
α
OC
(
1
−
N
sb,old
)
Δ
ε
ˉ
p
ε
0
,
OC
N 
sb,new
​
 =N 
sb,old
​
 +α 
OC
​
 (1−N 
sb,old
​
 ) 
ε 
0,OC
​
 
Δ 
ε
ˉ
  
p
 
​
 
The Olson–Cohen martensite fraction is
f
M
OC
=
1
−
exp
⁡
[
−
β
OC
N
sb,new
n
OC
]
f 
M
OC
​
 =1−exp[−β 
OC
​
 N 
sb,new
n 
OC
​
 
​
 ]
The JMAK fraction is
f
M
JMAK
=
1
−
exp
⁡
[
−
K
J
(
ε
ˉ
new
p
)
n
J
]
f 
M
JMAK
​
 =1−exp[−K 
J
​
 ( 
ε
ˉ
  
new
p
​
 ) 
n 
J
​
 
 ]
The combined base fraction is
f
M
base
=
1
−
(
1
−
f
M
OC
)
(
1
−
f
M
JMAK
)
f 
M
base
​
 =1−(1−f 
M
OC
​
 )(1−f 
M
JMAK
​
 )
Twin coupling is then added:
f
M
base
←
f
M
base
+
α
T
M
Δ
f
tw
(
1
−
f
M
base
)
f 
M
base
​
 ←f 
M
base
​
 +α 
TM
​
 Δf 
tw
​
 (1−f 
M
base
​
 )
The martensite fraction is enforced to be monotonic and bounded:
f
M
old
≤
f
M
new
≤
1
f 
M
old
​
 ≤f 
M
new
​
 ≤1
The increment is also capped for explicit stability.
Step 11: Re-check twin fraction
The twin fraction is reduced if necessary to satisfy
f
M
new
+
f
tw
new
≤
1
f 
M
new
​
 +f 
tw
new
​
 ≤1
Step 12: Assemble total inelastic strain
Δ
ε
inel
=
Δ
ε
cp
+
Δ
ε
tw
Δε 
inel
 =Δε 
cp
 +Δε 
tw
 
Step 13: Stress update
The elastic strain increment is
Δ
ε
e
=
Δ
ε
−
Δ
ε
inel
Δε 
e
 =Δε−Δε 
inel
 
The stress increment is
Δ
σ
=
D
e
:
Δ
ε
e
Δσ=D 
e
 :Δε 
e
 
Then the martensite volumetric eigenstrain correction is applied:
Δ
σ
i
i
←
Δ
σ
i
i
−
K
Δ
f
M
Δ
V
/
V
Δσ 
ii
​
 ←Δσ 
ii
​
 −KΔf 
M
​
 ΔV/V
Step 14: State variable update
All internal variables are written to stateNew.
21. Material Property Array
The subroutine expects:
1
Index
Symbol
Meaning
Typical unit
1
E
Young’s modulus
stress unit, e.g. MPa
2
nu
Poisson’s ratio
—
3
gam0
Reference slip rate
1/s
4
mexp
Rate-sensitivity exponent
—
5
tau0
Initial CRSS
stress unit
6
taus
Saturation CRSS
stress unit
7
h0
Voce hardening modulus
stress unit
8
alpha_OC
Olson–Cohen shear-band generation coefficient
—
9
beta_OC
Olson–Cohen nucleation coefficient
—
10
n_OC
Olson–Cohen exponent
—
11
eps0_OC
Reference plastic strain for Olson–Cohen
—
12
K_J
JMAK rate constant
—
13
n_J
JMAK/Avrami exponent
—
14
T_ref
Reference temperature
K
15
dVoV
Volumetric transformation strain 
Δ
V
/
V
ΔV/V
—
16
tau_twin
Critical resolved shear stress for twinning
stress unit
17
gam_twin
Characteristic twin shear
—
18
A_tw
Maximum twin volume fraction
—
19
B_tw
TANH transition width
stress unit
20
alpha_TM
Twin-to-martensite coupling coefficient
—
22. State Variables
The subroutine requires:
1
State variable
Index
Meaning
1–12
sv(1:12)
Accumulated slip shear strain on slip systems 1–12
13–24
sv(13:24)
Current CRSS on slip systems 1–12
25
sv(25)
Martensite volume fraction 
f
M
f 
M
​
 
26
sv(26)
Twin volume fraction 
f
tw
f 
tw
​
 
27
sv(27)
Equivalent plastic strain from slip
28
sv(28)
Normalized shear-band density 
N
sb
N 
sb
​
 
29
sv(29)
Accumulated transformation strain
30
sv(30)
Temperature
These can be requested as Abaqus field outputs using SDV1 through SDV30.
Recommended outputs:
Output
SDV
Use
Martensite fraction
25
Phase evolution
Twin fraction
26
TWIP evolution
Equivalent plastic strain
27
Plastic localization
Shear-band density
28
Transformation nucleation
Transformation strain
29
Volumetric transformation effect
Slip-system CRSS
13–24
Hardening evolution
23. Required Abaqus Input Deck
A minimal material definition should contain:
1234567891011
Example placeholder:
1234567891011
The above values are placeholders only and must be calibrated for the specific steel grade.
24. Crystallographic Orientation
The slip and twin systems are defined in the local material coordinate system. Therefore, for polycrystalline simulations, each grain or grain group must be assigned an orientation.
Example:
12
For Euler-angle-based orientation assignment, use the appropriate Abaqus orientation definition or preprocessing tool.
For a Voronoi RVE:
create one element set per grain,
assign one orientation per grain,
assign the same TRIP VUMAT material to all grains,
let Abaqus rotate the material frame grain by grain.
The VUMAT itself does not store or rotate grain orientations internally.
25. Voigt Component Ordering
Abaqus/Explicit VUMAT uses tensorial shear strain components, not engineering shear strains.
The helper routines in this subroutine assume the following order:
1
This is different from some Abaqus/Standard UMAT conventions. Users modifying the tensor conversion routines should be careful.
26. Numerical Stability Features
Explicit simulations can become unstable if internal variables evolve too rapidly. This VUMAT includes several protective features.
26.1 Power-law cap
The power-law term is capped:
∣
τ
g
∣
1
/
m
≤
10
8
​
  
g
τ
​
  
​
  
1/m
 ≤10 
8
 
This prevents overflow for very high stress ratios.
26.2 Slip increment cap
Each slip increment is limited:
∣
Δ
γ
α
∣
≤
0.01
∣Δγ 
α
 ∣≤0.01
This prevents excessively large slip in one increment.
26.3 Anti-overshoot limiter
The total deviatoric slip strain increment cannot exceed the imposed deviatoric strain increment. This avoids stress reversal due to over-correction.
26.4 Twin fraction rate cap
Δ
f
tw
≤
0.02
Δf 
tw
​
 ≤0.02
26.5 Martensite fraction rate cap
Δ
f
M
≤
0.02
Δf 
M
​
 ≤0.02
26.6 Argument clamping
Exponential and hyperbolic tangent arguments are clamped to safe ranges to prevent overflow or numerical singularities.
26.7 Monotonic phase evolution
Martensite and twin fractions are prevented from decreasing during unloading in the current implementation.
27. Strengths of the Implemented Model
27.1 Mechanism-based description
The model directly represents the main deformation and transformation mechanisms in TRIP/TWIP steels:
slip,
twinning,
martensitic transformation,
transformation dilatation,
twin-assisted transformation.
27.2 Crystallographic resolution
Because slip and twinning are defined on crystallographic systems, the model can represent orientation-dependent behavior and can be used in polycrystalline RVE simulations.
27.3 Robust explicit implementation
The subroutine includes several stabilization features that are essential for Abaqus/Explicit:
slip caps,
volume-fraction caps,
anti-overshoot limiter,
smooth TANH kinetics,
bounded phase fractions.
27.4 Flexible transformation kinetics
The combination of Olson–Cohen and JMAK kinetics gives flexibility in fitting experimental phase fraction evolution.
27.5 Twin–martensite interaction
The coupling term allows the model to capture the experimentally observed effect that twins may promote martensite nucleation.
27.6 Practical calibratability
The model uses a manageable number of physically meaningful parameters. It can be calibrated using:
uniaxial tension,
XRD phase fraction measurements,
EBSD texture data,
dilatometry,
digital image correlation,
strain-rate sensitivity tests.
28. Limitations and Modeling Simplifications
The model is intentionally designed as a robust engineering model. Users should understand its simplifications.
28.1 Isotropic elasticity
Elastic anisotropy of FCC austenite is not included. This is acceptable for many engineering simulations but may affect texture-sensitive elastic response.
28.2 Self-hardening only
Latent hardening between slip systems is not currently included. This can affect multislip hardening and texture evolution.
28.3 Martensite is not a separate crystal plasticity phase
Martensite is represented by volume fraction and eigenstrain, not by its own BCT/BCC slip systems. Therefore, the model does not directly capture martensite plasticity or two-phase load partitioning at the crystal level.
28.4 No explicit variant selection
The current transformation eigenstrain is isotropic volumetric expansion. Martensite variant selection and transformation texture are not modeled.
28.5 No explicit stress triaxiality or Lode angle dependence
The current kinetic laws do not explicitly include stress triaxiality or Lode angle. These can be added in future versions.
28.6 Temperature is fixed by default
The model stores temperature but does not update it from plastic dissipation. Adiabatic heating is not automatically included.
28.7 Twinning is irreversible and unidirectional
Detwinning, twin reorientation, and predominant twin reorientation are not included.
28.8 Energy output is approximate
The internal energy and inelastic energy outputs are not fully thermodynamically consistent. They should not be used for precise energy balance studies unless the energy update is improved.
28.9 Forward explicit update
The constitutive update is not a fully implicit return-mapping algorithm. It is designed for Abaqus/Explicit with small stable increments.
29. Comparison with Macroscopic Constitutive Models
Feature
Macroscopic 
J
2
J 
2
​
 /Hill/Barlat model
Present crystal plasticity TRIP model
Crystallographic slip
No
Yes
Texture sensitivity
Limited or empirical
Yes
Slip-system hardening
No
Yes
Twin systems
No
Yes
Phase fraction evolution
Usually no
Yes
Transformation eigenstrain
Usually no
Yes
Twin–martensite coupling
No
Yes
Computational cost
Low
Moderate to high
Calibration complexity
Lower
Higher
Physical resolution
Low
High
30. Recommended Calibration Strategy
30.1 Elastic parameters
Obtain from literature or experiments:
E
E,
ν
ν.
For steel, typical values are approximately:
E
≈
190
–
210
 GPa
E≈190–210 GPa
ν
≈
0.28
–
0.30
ν≈0.28–0.30
Use consistent Abaqus units.
30.2 Slip parameters
Calibrate using single-crystal or polycrystalline stress–strain data:
tau0 from initial yield,
taus from saturation stress,
h0 from initial hardening slope,
gam0 and mexp from strain-rate sensitivity.
For a rough initial estimate:
σ
yield
≈
τ
0
M
Schmid,max
σ 
yield
​
 ≈ 
M 
Schmid,max
​
 
τ 
0
​
 
​
 
where 
M
Schmid,max
M 
Schmid,max
​
  is the maximum Schmid factor for the loading direction.
30.3 Transformation parameters
Use experimental martensite fraction versus plastic strain data.
Fit:
alpha_OC,
beta_OC,
n_OC,
eps0_OC,
K_J,
n_J.
The Olson–Cohen parameters mainly control the early transformation and shear-band nucleation behavior. The JMAK parameters control the overall sigmoidal evolution.
30.4 Volumetric transformation strain
Use crystallographic data or dilatometry.
Typical values:
Δ
V
/
V
≈
0.02
–
0.04
ΔV/V≈0.02–0.04
30.5 Twinning parameters
Calibrate using twin fraction measurements or mechanical response:
tau_twin from twinning onset,
gam_twin from crystallographic twin shear, usually about 0.707 for FCC,
A_tw from maximum observed twin fraction,
B_tw from the sharpness of twin activation.
30.6 Twin–martensite coupling
Set alpha_TM = 0 initially. Then increase it if experimental evidence shows that twinning accelerates martensite formation.
31. Suggested Validation Cases
31.1 Single-element uniaxial tension
Check:
elastic slope,
initial yield,
hardening shape,
stability.
31.2 Single-crystal orientation test
Apply tension along known crystallographic directions and compare yield stress against Schmid-factor predictions.
31.3 Polycrystalline RVE
Assign random or measured orientations and examine:
macroscopic stress–strain response,
phase fraction evolution,
strain localization,
grain-level heterogeneity.
31.4 Phase fraction evolution
Compare predicted 
f
M
f 
M
​
  versus strain against:
XRD,
EBSD,
magnetic measurements,
neutron diffraction.
31.5 Twin fraction evolution
Compare predicted 
f
tw
f 
tw
​
  against:
EBSD twin boundary analysis,
TEM observations,
optical microscopy.
31.6 Strain-rate sensitivity
Run simulations at different strain rates and verify that the power-law slip rule produces realistic rate dependence.
32. Possible Extensions
The model can be extended in several directions.
32.1 Latent hardening
Introduce a hardening matrix:
g
˙
α
=
∑
β
h
α
β
∣
γ
˙
β
∣
g
˙
​
  
α
 = 
β
∑
​
 h 
αβ
 ∣ 
γ
˙
​
  
β
 ∣
This improves multislip and texture predictions.
32.2 Cubic elastic anisotropy
Replace isotropic Hooke elasticity with cubic elasticity using:
C
11
,
C
12
,
C
44
C 
11
​
 ,C 
12
​
 ,C 
44
​
 
32.3 Full finite-strain hyperelastic formulation
Implement a multiplicative decomposition with an explicit update of:
F
e
,
F
p
,
F
tr
,
F
tw
F 
e
 ,F 
p
 ,F 
tr
 ,F 
tw
 
32.4 Stress-state-dependent transformation kinetics
Introduce triaxiality and Lode angle dependence:
K
J
=
K
J
(
η
,
θ
ˉ
)
K 
J
​
 =K 
J
​
 (η, 
θ
ˉ
 )
or
α
OC
=
α
OC
(
η
,
θ
ˉ
)
α 
OC
​
 =α 
OC
​
 (η, 
θ
ˉ
 )
32.5 Thermomechanical coupling
Update temperature from plastic work:
Δ
T
=
β
Taylor
ρ
c
p
σ
:
Δ
ε
p
ΔT= 
ρc 
p
​
 
β 
Taylor
​
 
​
 σ:Δε 
p
 
where 
β
Taylor
β 
Taylor
​
  is the Taylor–Quinney coefficient.
32.6 Two-phase crystal plasticity
Represent martensite as a separate phase with its own slip systems:
FCC austenite slip,
BCT/BCC martensite slip,
phase-specific hardening,
load partitioning.
32.7 Transformation variant selection
Introduce martensite variants and transformation shear tensors:
ε
tr
=
∑
v
f
v
ε
tr
,
v
ε 
tr
 = 
v
∑
​
 f 
v
 ε 
tr,v
 
32.8 Twin reorientation
Implement predominant twin reorientation, PTR, or crystallographic reorientation of twinned regions.
32.9 Damage and fracture
Couple the model with continuum damage mechanics to predict ductile fracture after extensive transformation and twinning.
33. Model Selection Summary
Model component
Selected model
Main reason
Elasticity
Isotropic Hooke
Efficient, stable, sufficient for small elastic strains
Plastic flow
Crystal plasticity slip systems
Captures crystallographic mechanisms directly
Slip systems
FCC {111}<110>
Standard austenite slip family
Slip kinetics
Power-law viscoplasticity
Rate dependence, explicit stability
Slip hardening
Voce per system
Captures saturation hardening
Macroscopic yield surface
Not used
Yielding emerges from slip-system CRSS
Martensite kinetics
Olson–Cohen + JMAK
Combines shear-band nucleation and sigmoidal growth
Transformation strain
Isotropic volumetric eigenstrain
Captures dominant volume expansion robustly
Twinning systems
FCC {111}<112>
Standard FCC twin family
Twinning criterion
Polar CRSS
Captures directional nature of twinning
Twinning kinetics
TANH law
Smooth, stable, saturating
Twin–martensite coupling
Linear coupling coefficient
Captures twin-assisted martensite nucleation
Numerical stabilization
Caps and limiters
Required for robust Abaqus/Explicit simulations
34. Practical Usage Recommendations
Use consistent units.
If stress is in MPa, density and geometry must be consistent with MPa.
Start with a single-element test.
Verify stability before using the model in a large RVE.
Use small stable time increments.
Explicit simulations are sensitive to rapid internal variable changes.
Monitor phase fractions.
Ensure 
f
M
f 
M
​
  and 
f
tw
f 
tw
​
  remain physically reasonable.
Check orientation assignment.
Incorrect orientations can produce unrealistic texture and yield behavior.
Calibrate transformation and twinning separately if possible.
Use phase fraction data for transformation and twin fraction data for twinning.
Disable coupling initially.
Set alpha_TM = 0 first, then enable it after the base model is calibrated.
35. Known Items to Verify Before Publication
Before using this subroutine for publication or industrial prediction, verify:
Initial yield stress against Schmid-factor prediction.
Hardening curve against experimental tensile data.
Martensite fraction evolution against XRD/EBSD data.
Twin fraction evolution against microstructure measurements.
Strain-rate sensitivity against split-Hopkinson or multi-rate tests.
Sensitivity to explicit time increment and mass scaling.
Orientation dependence in single-crystal benchmark tests.
Volume expansion effect under constrained boundary conditions.
Energy outputs if energy balance is important.
Twin-system sign convention relative to the intended crystallographic texture.
36. Summary
This VUMAT provides a robust, mechanism-based constitutive model for TRIP/TWIP steels in Abaqus/Explicit. It combines:
rate-dependent FCC crystal plasticity,
Voce hardening,
Olson–Cohen shear-band transformation nucleation,
JMAK strain-driven transformation kinetics,
isotropic transformation volumetric expansion,
polar CRSS twinning activation,
smooth TANH twin kinetics,
twin-assisted martensite nucleation,
explicit numerical stabilization.
The model is not intended to be a fully resolved two-phase micromechanical model. Instead, it is a computationally efficient, physically motivated material-point model capable of capturing the essential coupling between slip, twinning, and martensitic transformation in metastable austenitic/TRIP steels.
It is especially suitable for simulations where the simultaneous evolution of plasticity, phase transformation, and twinning controls the macroscopic mechanical response.
37. Key References for Further Reading
J. R. Rice, crystal plasticity foundations.
D. Peirce, R. J. Asaro, A. Needleman, rate-dependent crystal plasticity.
S. R. Kalidindi, C. A. Bronkhorst, L. Anand, crystal plasticity calibration for FCC metals.
Y. Huang, Harvard crystal plasticity VUMAT/UMAT implementation.
G. B. Olson, M. Cohen, strain-induced martensitic transformation.
D. P. Koistinen, R. E. Marburger, martensite transformation kinetics.
W. A. Johnson, R. F. Mehl, M. Avrami, A. N. Kolmogorov, JMAK transformation theory.
Y. Tomita, T. Iwamoto, TRIP steel constitutive modeling.
P. Van Houtte, predominant twin reorientation.
Advanced high-strength steel literature for TRIP/TWIP mechanisms.
38. License and Usage Note
This model is intended for research and educational use. Industrial application requires careful calibration and validation against experimental data for the specific steel grade, heat treatment, texture, strain rate, and temperature range.
