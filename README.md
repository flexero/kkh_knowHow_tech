# 01 · HVAC Design — Complete Implementation Reference

> Heating, Ventilation and Air-Conditioning design: principles, standard equations, real design
> parameters, a step-by-step method and a fully worked example.
> Companion files: `hvac_engine.py` (Python stdlib engine) · `index.html` (offline visual guide + calculators).
> **Run offline:** `run.bat` (Windows) / `sh run.sh` (Linux/macOS) or `python hvac_engine.py`. Python 3.6–3.15, standard library only: nothing to install (see the root README).
> Units: SI throughout. Main references: ASHRAE Handbook – Fundamentals (2017/2021), ASHRAE 55, 62.1, 90.1,
> SMACNA HVAC Duct Design, CIBSE Guide A/C. Always check against your **national code**
> (e.g. Vietnam: TCVN 5687:2010 ventilation & AC design, QCVN 09:2017/BXD energy efficiency).

---

**How this book is organised:** **Part A** gives the scientific fundamentals (why). **Part B** is the design
procedure (how to size the system). **Part C** covers equipment engineering and manufacturing (how the
hardware is designed, built and tested). **Part D** covers standards, QA and operation. Every numbered
example is reproduced by `python hvac_engine.py`.

## Contents
**Part A — Scientific fundamentals**
- [A1 Units & constants](#a1-units-constants-and-conversions) · [A2 Thermodynamics](#a2-thermodynamics) · [A3 Moist air](#a3-moist-air-psychrometric-theory) · [A4 Heat transfer](#a4-heat-transfer) · [A5 Fluid mechanics](#a5-fluid-mechanics) · [A6 Refrigeration cycles](#a6-refrigeration-cycles)

**Part B — Design procedure**

1. [Design workflow](#1-design-workflow)
2. [Indoor & outdoor design conditions](#2-indoor--outdoor-design-conditions)
3. [Psychrometrics](#3-psychrometrics)
4. [Heat transfer through the envelope](#4-heat-transfer-through-the-envelope)
5. [Solar radiation](#5-solar-radiation)
6. [Cooling load calculation](#6-cooling-load-calculation)
7. [Heating load calculation](#7-heating-load-calculation)
8. [Ventilation — ASHRAE 62.1](#8-ventilation--ashrae-621)
9. [Air-side process: supply air, mixing, coil](#9-air-side-process-supply-air-mixing-coil)
10. [System selection](#10-system-selection)
11. [Air distribution & diffusers](#11-air-distribution--diffusers)
12. [Duct design](#12-duct-design)
13. [Hydronic (pipe) design](#13-hydronic-pipe-design)
14. [Fans & pumps](#14-fans--pumps)
15. [Central plant: chillers, towers, boilers](#15-central-plant-chillers-towers-boilers)
16. [Energy, efficiency & controls](#16-energy-efficiency--controls)
17. [Worked example (full output)](#17-worked-example-full-output)
18. [Engine API reference](#18-engine-api-reference)
19. [Rules of thumb (for checking only)](#19-rules-of-thumb-for-checking-only)
20. [Design checklist](#20-design-checklist)
21. [Equipment selection & bill of materials (procurement)](#21-equipment-selection--bill-of-materials-procurement)
22. [Installation tools & test (TAB) instruments](#22-installation-tools--test-tab-instruments)

**Part C — Equipment engineering & manufacturing**
- [C1 Cooling coil design & fabrication](#c1-finned-tube-cooling-coil-design-rating-and-fabrication) · [C2 Fans](#c2-fans-engineering--construction) · [C3 Pumps](#c3-pumps-engineering--construction) · [C4 Duct construction & leakage](#c4-ductwork-construction--leakage) · [C5 AHU construction & filters](#c5-air-handling-unit-construction-en-1886--en-13053) · [C6 Insulation](#c6-thermal-insulation-engineering) · [C7 Acoustics](#c7-acoustics) · [C8 Vibration isolation](#c8-vibration-isolation) · [C9 DX/VRF installation](#c9-dx--vrf-refrigerant-installation-practice)

**Part D — Standards, QA, operation**
- [D1 Standards map](#d1-standards-map) · [D2 Inspection & test plan](#d2-inspection--test-plan-itp-with-acceptance-criteria) · [D3 Operation & maintenance](#d3-operation--maintenance-preventive-schedule)

---

# PART A — SCIENTIFIC FUNDAMENTALS

Part A gives the physics behind every design rule in Parts B–D. Each formula used later is derived or
stated here, with its validity range.

## A1. Units, constants and conversions

| Quantity | Value |
|---|---|
| Universal gas constant R̄ | 8.314462 J/(mol·K) |
| Dry air: M = 28.966 g/mol, R_da | 287.042 J/(kg·K) |
| Water vapour: M = 18.015 g/mol, R_w | 461.524 J/(kg·K) |
| c_p dry air / water vapour / liquid water | 1.006 / 1.86 / 4.186 kJ/(kg·K) |
| h_fg water at 0 °C / 25 °C / 100 °C | 2501 / 2442 / 2257 kJ/kg |
| Stefan–Boltzmann σ | 5.670374 × 10⁻⁸ W/(m²·K⁴) |
| Standard gravity g | 9.80665 m/s² |
| Standard atmosphere | 101.325 kPa |
| 1 TR (ton of refrigeration) | 3.517 kW = 12 000 Btu/h |
| 1 kW | 3412 Btu/h |
| 1 m³/s | 1000 L/s = 3600 m³/h = 2119 CFM |
| 1 Pa | 0.004015 in.wg; 1 in.wg = 249.1 Pa |
| 1 bar | 100 kPa = 10.2 m water column |
| 1 W/(m²·K) | 0.1761 Btu/(h·ft²·°F) |
| 1 m²·K/W (R-value) | 5.678 h·ft²·°F/Btu |

## A2. Thermodynamics

**First law, steady-flow open system** (every HVAC component: coil, fan, pump, chiller, mixing box):
```
Q̇ − Ẇ = Σ ṁ_out (h + V²/2 + g z) − Σ ṁ_in (h + V²/2 + g z)
```
In HVAC the kinetic and potential terms are usually negligible, so:
- Coil: `Q̇ = ṁ_a (h_out − h_in) + ṁ_condensate · h_f`. The condensate term is about 0.5 % and is normally dropped.
- Fan: `Ẇ = ṁ Δp_t / (ρ η)`; all fan power ends up as heat in the air stream.
- Adiabatic mixing: `Σṁ h` is conserved, and so is `Σṁ W` (mass of water).

**Second law & refrigeration limits.** No cycle can beat the reversible Carnot cycle between T_L and T_H (K):
```
COP_cooling,Carnot = T_L / (T_H − T_L)        COP_heating,Carnot = T_H / (T_H − T_L)
Second-law efficiency  η_II = COP_real / COP_Carnot      (good chillers 0.55–0.70)
```
Example: evaporating 0 °C (273.15 K), condensing 40 °C (313.15 K) gives COP_Carnot = 273.15/40 = **6.83**.
A real cycle at COP 4.5 has η_II = 0.66.
**Each 1 K rise in evaporating temperature or 1 K fall in condensing temperature improves COP by about 2–4 %.**
This is the physics behind chilled-water reset and condenser-water reset.

**Exergy** (available work) of heat at T relative to ambient T₀: `Ex = Q (1 − T₀/T)`. This explains why
low-temperature-lift systems (radiant cooling at 16 °C, heat pumps at 35 °C) are efficient.

## A3. Moist air (psychrometric theory)

Moist air is treated as a mixture of two ideal gases (error < 0.7 % below 50 °C).
- **Dalton's law:** `p = p_a + p_w`.
- **Humidity ratio, derived:** `W = m_w/m_a = (p_w V/R_w T)/(p_a V/R_a T) = (R_a/R_w)·p_w/(p − p_w)`, with
  `R_a/R_w = 287.042/461.524 = 0.621945`.
- **Enthalpy per kg of dry air:** `h = c_pa·t + W(h_fg,0 + c_pw·t) = 1.006 t + W(2501 + 1.86 t)`.
  The reference state is h = 0 for dry air at 0 °C and liquid water at 0 °C.
- **Adiabatic saturation (thermodynamic wet bulb t\*)**, energy balance on a long adiabatic channel:
  `h + (W_s* − W)·h_w* = h_s*`. Solving for W gives the wet-bulb equation used in the engine.
  A psychrometer wet bulb equals t\* within about 0.3 K when the air speed is ≥ 3 m/s (the Lewis number of air ≈ 1).

**Origin of the practical constants (sea level):**
```
Sensible:  Q = ρ c_p V Δt = 1.20 kg/m³ × 1.026 kJ/kgK × V[m³/s] × Δt  → 1.23 W per (L/s·K)
Latent:    Q = ρ h_fg V ΔW = 1.20 × 2501 × V × ΔW                    → 3010 W per (L/s · kg/kg)
Total:     Q = ρ V Δh                                                 → 1.20 W per (L/s · kJ/kg)
```
c_p = 1.026 includes the vapour (1.006 + 1.86 × 0.010). At altitude, ρ ∝ p, so multiply by p/101.325.

## A4. Heat transfer

### A4.1 Conduction (Fourier: q″ = −k dT/dx)

| Geometry | Thermal resistance R [K/W] |
|---|---|
| Plane wall, thickness L, area A | `L / (k A)` |
| Cylinder shell r₁→r₂, length ℓ | `ln(r₂/r₁) / (2π k ℓ)` |
| Sphere shell | `(1/r₁ − 1/r₂) / (4π k)` |
| Convection surface | `1 / (h A)` |
| Contact / fouling | `R″_f / A` |

Resistances in series add; parallel paths add as conductances (1/R). **Critical radius of insulation** on a
cylinder is `r_c = k/h`. Insulating a pipe smaller than r_c *increases* heat flow. For k = 0.036 and
h = 9, r_c = 4 mm, so this is irrelevant for normal pipes but matters for thin wires.

**Transient conduction / thermal mass.** For a lumped body (Biot number `Bi = hL_c/k < 0.1`):
`T − T∞ = (T₀ − T∞)·exp(−t/τ)`, with `τ = ρ c V / (h A)`. This is the physics behind the decrement factor,
time lag and radiant time series in §4.4 and §6. Heavy walls have long τ, so heat peaks are delayed and damped.
For a periodic (24 h) surface temperature, a semi-infinite slab damps the amplitude as
`exp(−x/δ)` with penetration depth `δ = √(2α/ω)`. Here α = k/(ρc) is the thermal diffusivity and ω = 2π/86400 s⁻¹.
For concrete (α ≈ 7 × 10⁻⁷ m²/s), δ ≈ 0.14 m: a 0.3 m concrete wall damps the swing to about e^(−2.1) ≈ 0.12
and delays it by about x/δ radians ≈ 8 h. This matches f ≈ 0.15–0.25 and φ ≈ 9–10 h in the §4.4 table.

### A4.2 Convection (Newton: q = h A ΔT)

| Situation | h [W/(m²·K)] |
|---|---|
| Free convection, air | 2–25 (indoor surfaces 3–8) |
| Forced convection, air | 25–250 (finned coil 40–80) |
| Forced convection, water | 500–10 000 (coil tubes 2000–6000) |
| Condensing refrigerant | 2000–25 000 |
| Boiling (evaporating) refrigerant | 2500–35 000 |

**Correlations** (Nu = hD/k, Re = ρVD/μ, Pr = μc_p/k):
- Tube, laminar, fully developed, constant wall temperature: `Nu = 3.66`.
- Tube, turbulent, **Dittus–Boelter**: `Nu = 0.023 Re^0.8 Pr^n` (n = 0.4 heating the fluid, 0.3 cooling), valid for Re > 10⁴.
- Tube, turbulent, **Gnielinski** (more accurate; 3000 < Re < 5×10⁶; used in the engine):
  `Nu = (f/8)(Re − 1000)Pr / [1 + 12.7 (f/8)^0.5 (Pr^(2/3) − 1)]`, with `f = (0.790 ln Re − 1.64)⁻²`.
- Flat plate, laminar: `Nu_L = 0.664 Re_L^0.5 Pr^(1/3)`; turbulent: `0.037 Re_L^0.8 Pr^(1/3)`.
- Vertical plate, natural convection: `Nu = 0.59 Ra^(1/4)` for 10⁴ < Ra < 10⁹, and `0.10 Ra^(1/3)` above;
  `Ra = gβΔT L³/(να)`.
- Finned-tube coil, air side: Gray & Webb, see §C1.

### A4.3 Radiation

`q = ε σ A (T₁⁴ − T₂⁴)` for a small grey surface in a large enclosure. Linearised: `h_r = 4 ε σ T_m³`
≈ 5.1 W/(m²·K) at 20 °C with ε = 0.9. Interior surface combined coefficient: h_c + h_r ≈ 7.7 W/(m²·K),
which gives R_si = 0.13. Emissivity: building materials and paint 0.85–0.95; aluminium foil 0.03–0.10;
low-e glass coating 0.02–0.15.

### A4.4 Fins

Straight fin, adiabatic tip: `η_f = tanh(mL)/(mL)`, with `m = √(hP/(k A_c)) ≈ √(2h/(k t))` for thin fins.
Surface efficiency of the finned side: `η_o = 1 − (A_fin/A_total)(1 − η_f)`.
Plate fins on tubes use the Schmidt equivalent circular fin (§C1). Wet fins conduct heat *and* condense
water, so the effective h is h·c_s/c_p (c_s = dh_sat/dT), which lowers fin efficiency (Threlkeld).

### A4.5 Heat exchangers

```
Overall: 1/(UA) = 1/(h_i A_i) + R″_fi/A_i + R_wall + R″_fo/A_o + 1/(η_o h_o A_o)
LMTD method:  Q = U A F ΔT_lm,   ΔT_lm = (ΔT₁ − ΔT₂)/ln(ΔT₁/ΔT₂)
ε-NTU method: C = ṁ c_p, C_r = C_min/C_max, NTU = UA/C_min, Q = ε C_min (T_h,in − T_c,in)
  counter-flow : ε = [1 − e^(−NTU(1−C_r))] / [1 − C_r e^(−NTU(1−C_r))]    (C_r = 1: ε = NTU/(1+NTU))
  parallel     : ε = [1 − e^(−NTU(1+C_r))] / (1 + C_r)
  cross-flow, both unmixed : ε ≈ 1 − exp{ (1/C_r) NTU^0.22 [exp(−C_r NTU^0.78) − 1] }
  C_r = 0 (condenser/evaporator) : ε = 1 − e^(−NTU)
```
Coils with ≥ 4 rows and counter-flow circuiting behave close to pure counter-flow.

| Fouling factor R″_f | m²·K/W |
|---|---|
| Closed-loop treated chilled water | 0.000044–0.000088 |
| Cooling-tower water, treated | 0.000176 |
| City water | 0.000176 |
| River / untreated | 0.00035–0.00053 |
| Refrigerant | 0.000018 |

AHRI 550/590 rates chillers with evaporator fouling 0.000018 m²·K/W and condenser fouling 0.000044 m²·K/W.

## A5. Fluid mechanics

- **Continuity:** `ṁ = ρ A V = constant`.
- **Energy with losses** (Bernoulli extended), in pressure units:
  `p₁ + ½ρV₁² + ρgz₁ = p₂ + ½ρV₂² + ρgz₂ + Δp_loss`.
  In ducts, total pressure = static + velocity pressure, and `p_v = ½ρV² ≈ 0.602 V²` [Pa].
- **Reynolds number** `Re = VD/ν`: laminar below 2300, turbulent above about 4000.
  HVAC ducts run at Re ≈ 10⁵ and water pipes at 10⁴–10⁶, so they are always turbulent.
- **Darcy–Weisbach** `Δp = f (L/D) ½ρV²`. Laminar flow: f = 64/Re. Turbulent flow: Colebrook (implicit),
  or the explicit Swamee–Jain form `f = 0.25/[log₁₀(ε/3.7D + 5.74/Re^0.9)]²` (within ±1 %).
- **Minor losses:** `Δp = K ½ρV²`, or as an equivalent length `L_e = K D / f`.
- **Hydraulic diameter** of a non-circular duct: `D_h = 4A/P`. Use D_h for Re; for equal friction *and* equal
  flow use the Huebscher D_e (§12.2).
- **System curve:** `Δp = Δp_static + k V̇²`. The operating point is where it intersects the fan or pump curve.
- **Euler turbomachine equation:** `Δp_t = ρ (u₂ c_u2 − u₁ c_u1)`, or for pumps `H = (u₂ c_u2 − u₁ c_u1)/g`,
  where u is blade speed and c_u is the tangential absolute velocity.
- **Dimensional analysis** gives the affinity laws `V̇ ∝ n D³`, `Δp ∝ ρ n² D²`, `P ∝ ρ n³ D⁵`.
- **Specific speed** `n_q = n √V̇ / H^0.75` (rpm, m³/s, m). Radial impellers: 10–40; mixed flow: 40–160;
  axial: > 160. Use it to choose the impeller type.
- **Water hammer** (Joukowsky): `Δp = ρ a ΔV`, with wave speed a ≈ 1200–1400 m/s in steel pipe.
  Stopping 2 m/s instantly gives about 26 bar. Valves must close slower than `2L/a` and pumps need soft stop.

## A6. Refrigeration cycles

### A6.1 Vapour-compression cycle (p–h diagram)

```
 p ▲         3 ●───────────────● 2s/2      1→2  compression (η_is = (h2s−h1)/(h2−h1) = 0.65–0.80)
   │  subcool  │   condenser    ╱          2→3  condensation + subcooling (rejects q_c = h2 − h3)
   │           │               ╱           3→4  throttling, h4 = h3 (isenthalpic)
   │         4 ●───────────────● 1          4→1  evaporation + superheat (refrigerating effect q_e = h1 − h4)
   │               evaporator   superheat
   └──────────────────────────────────────▶ h
 ṁ = Q_e / (h1 − h4)    W = ṁ (h2 − h1)    COP = (h1 − h4)/(h2 − h1)    Q_c = Q_e + W
 Displacement V̇_d = ṁ · v1 / η_vol     (η_vol ≈ 0.85–0.95 scroll/screw, falls with pressure ratio)
```

**Example: R134a, T_e = 0 °C with 5 K superheat, T_c = 40 °C with 5 K subcooling, Q_e = 100 kW.**
Table values: h1 = 402.7, h3 = h4 = 249.7, h2s = 427.0 kJ/kg, v1 = 0.0709 m³/kg (IIR reference, approximate;
verify with NIST REFPROP or CoolProp). With η_is = 0.72, `vapour_compression()` gives:
h2 = 436.5 kJ/kg, q_e = 153.0 kJ/kg, ṁ = 0.654 kg/s, W = 22.1 kW, Q_c = 122.1 kW, **COP = 4.53**
(η_II = 4.53/6.83 = 0.66), displacement = 185 m³/h.

### A6.2 Refrigerants

| Refrigerant | Type | GWP (AR4) | Safety class | NBP °C | T_crit °C | Typical use |
|---|---|---|---|---|---|---|
| R32 | HFC | 675 | A2L | −51.7 | 78.1 | split, VRF, small chillers |
| R410A | HFC blend | 2088 | A1 | −51.4 | 71.3 | split/VRF (being phased down) |
| R134a | HFC | 1430 | A1 | −26.1 | 101.1 | screw/centrifugal chillers |
| R513A | HFO blend | 631 | A1 | −29.2 | ≈ 95 | R134a replacement |
| R1234ze(E) | HFO | < 1 | A2L | −19.0 | 109.4 | low-GWP chillers |
| R1234yf | HFO | < 1 | A2L | −29.5 | 94.7 | automotive, chillers |
| R290 (propane) | HC | 3 | A3 | −42.1 | 96.7 | small units, heat pumps (charge limits!) |
| R717 (ammonia) | natural | 0 | B2L | −33.3 | 132.3 | industrial, large chillers |
| R744 (CO₂) | natural | 1 | A1 | −78.4 (subl.) | 31.0 | transcritical heat pumps, refrigeration |

Safety class (ISO 817 / ASHRAE 34): A = lower toxicity, B = higher toxicity; 1 = no flame propagation,
2L = mildly flammable, 2 = flammable, 3 = highly flammable. Regulation: Kigali Amendment HFC phase-down
(national schedules), EU F-gas.

### A6.3 Design temperatures & components

| Parameter | Typical value |
|---|---|
| Evaporating temperature, chiller (CHW 6 °C leaving) | 2–5 °C (approach 1–3 K) |
| Evaporating temperature, DX air coil | 5–9 °C |
| Condensing, water-cooled (CW 37 °C leaving) | 38–40 °C (approach 1–3 K) |
| Condensing, air-cooled | ambient + 12–18 K (45–55 °C) |
| Superheat at evaporator outlet | 4–8 K (TXV/EEV setpoint) |
| Subcooling at condenser outlet | 3–8 K |
| Compressor discharge temperature limit | ≤ 110–120 °C (oil breakdown) |

| Compressor | Range | Capacity control | Notes |
|---|---|---|---|
| Rotary | 1–10 kW | inverter | room splits |
| Scroll | 3–60 kW (tandem/trio to 200 kW) | inverter, staging, digital | few moving parts |
| Reciprocating | 1–200 kW | cylinder unloading | legacy |
| Screw (twin) | 100–2000 kW | slide valve or VSD | robust, oil injected |
| Centrifugal | 300–20 000 kW | inlet guide vanes, VSD; magnetic bearings (oil-free) | best full- and part-load COP; watch surge |

Expansion devices: capillary (fixed, small units), TXV (thermostatic, constant superheat), EEV (electronic;
best for inverter/VRF), float valves (flooded chillers). Evaporators: DX (refrigerant in tubes),
flooded shell-and-tube, falling film, brazed plate. Condensers: air-cooled finned coil, shell-and-tube
water-cooled, brazed plate, evaporative.

**Refrigerant piping design criteria (ASHRAE Refrigeration Handbook):**

| Line | Velocity | Pressure-drop limit |
|---|---|---|
| Suction, horizontal | ≥ 3.8 m/s (oil return) | ≈ 1 K saturated-temperature equivalent (max 2 K) |
| Suction riser | ≥ 7.6 m/s (oil return); double riser or traps if capacity is reduced | — |
| Discharge | 10–18 m/s (≥ 7.6 m/s in risers) | ≈ 1 K |
| Liquid | 0.5–1.5 m/s (avoid hammer and flashing) | ≤ 0.5–1 K; subcooling must exceed riser static loss |

### A6.4 Absorption & heat pumps

- **Single-effect LiBr–water absorption chiller:** COP ≈ 0.7, driven by 80–120 °C hot water or low-pressure steam.
- **Double-effect:** COP ≈ 1.2–1.4, driven by 150–180 °C steam or direct gas firing.
- **Constraints:** LiBr crystallisation limits the cooling-water temperature, and the cooling tower rejects
  Q_e(1 + 1/COP), about 2.4× Q_e versus 1.2× for electric chillers.
- **Heat pumps** (the same cycle used for heating): `COP_h = COP_c + 1`.

---

# PART B — DESIGN PROCEDURE

## 1. Design workflow

```
 ┌─────────────┐  ┌──────────────┐  ┌─────────────┐  ┌──────────────┐  ┌──────────────┐
 │1 Brief &    │→ │2 Design      │→ │3 Envelope & │→ │4 Zone loads  │→ │5 Ventilation │
 │  zoning     │  │  conditions  │  │  U-values   │  │  (hourly)    │  │  (62.1)      │
 └─────────────┘  └──────────────┘  └─────────────┘  └──────────────┘  └──────┬───────┘
                                                                              ↓
 ┌─────────────┐  ┌──────────────┐  ┌─────────────┐  ┌──────────────┐  ┌──────────────┐
 │10 Controls &│← │9 Plant       │← │8 Pipes &    │← │7 Ducts &     │← │6 Psychro:    │
 │  commission │  │  (chiller…)  │  │  pumps      │  │  fans        │  │  airflow,coil│
 └─────────────┘  └──────────────┘  └─────────────┘  └──────────────┘  └──────────────┘
```

| Step | Input | Output | Engine function |
|---|---|---|---|
| 1 Zoning | Architectural plans, use, schedules | Thermal zones (same orientation/use/schedule) | — |
| 2 Conditions | Climate data, comfort standard | Outdoor DB/MCWB/daily range, indoor T/RH | `AirState`, `outdoor_profile` |
| 3 Envelope | Wall/roof/glass build-ups | U, SHGC, decrement f, lag φ | `Opaque`, `Window` |
| 4 Loads | Zones + gains | Hourly sensible/latent, peak hour | `cooling_load`, `heating_load` |
| 5 Ventilation | Occupancy, area, space type | Outdoor air L/s | `outdoor_air`, `system_oa_multizone` |
| 6 Air process | Loads, supply ΔT, OA | Supply airflow, coil load, ADP, BF | `supply_airflow`, `mix`, `coil_adp_bf`, `coil_load` |
| 7 Ducts | Layout, flows | Sizes, ESP, fan selection | `size_duct_network`, `fan_power` |
| 8 Pipes | Coil loads, ΔT | Pipe DN, pump head | `water_flow`, `size_pipe`, `pump_power` |
| 9 Plant | Block load (coincident) | Chiller/tower/boiler size | `plant_sizing` |
| 10 Controls | Sequence of operation | Setpoints, resets, BMS points | §16 |

**Key principle — block vs. sum of peaks:** each zone/terminal is sized for its *own* peak hour;
AHUs serving several zones and the central plant are sized for the *coincident (block)* peak
(the hour when the sum is largest). Block load is typically 75–90 % of the sum of zone peaks.

---

## 2. Indoor & outdoor design conditions

### 2.1 Indoor (comfort) — ASHRAE 55 / ISO 7730

| Space | Summer DB °C | Summer RH % | Winter DB °C | Notes |
|---|---|---|---|---|
| Office, meeting | 23–26 (design 24) | 40–60 (design 50) | 20–22 | PMV −0.5…+0.5, PPD ≤ 10 % |
| Hotel room | 23–25 | 45–60 | 21–23 | |
| Retail | 24–26 | 45–60 | 19–21 | |
| Hospital ward | 24 | 40–60 | 22–24 | ASHRAE 170 |
| Operating room | 20–24 | 20–60 | 20–24 | ASHRAE 170, ≥20 ACH, positive |
| Data centre (inlet) | 18–27 | DP −9…15 °C | — | ASHRAE TC 9.9 recommended |
| Cleanroom (ISO 7) | 20–22 ±1 | 45 ±5 | — | Process dependent |

Air velocity in occupied zone: ≤ 0.20 m/s cooling (draught), ≤ 0.15 m/s heating.
Max humidity ratio for comfort (ASHRAE 55): 12 g/kg.

**PMV (Fanger) — operative temperature** `t_op ≈ (t_a + t_mrt)/2` when air speed < 0.2 m/s.
Clothing: summer 0.5 clo, winter 1.0 clo. Metabolic rate: office 1.1–1.2 met (1 met = 58.2 W/m²).

### 2.2 Outdoor design data

Use ASHRAE Climatic Design Conditions (or national standard):

| Purpose | Use |
|---|---|
| Cooling, comfort | 0.4 % or 1 % DB with **MCWB** (mean coincident wet bulb) + daily range |
| Dehumidification / OA units | 0.4 % or 1 % **DP** with MCDB (humid climates control latent load!) |
| Cooling tower | 0.4 % WB (evaporation design) |
| Heating | 99.6 % or 99 % DB |

`0.4 %` means the value is exceeded 0.4 % of the 8760 h year (≈ 35 h).

Example values used in this README (hot-humid, lat 21 °N, July): **DB 36.0 °C, MCWB 28.0 °C,
daily range 7 K**, winter 8 °C. Replace with your site's data.

### 2.3 Hourly outdoor temperature profile

`T(h) = T_DB,design − f(h) · DR`  (ASHRAE fraction of daily range, solar time):

| h | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| f | .87 | .92 | .96 | .99 | 1.00 | .98 | .93 | .84 | .71 | .56 | .39 | .23 |
| **h** | **13** | **14** | **15** | **16** | **17** | **18** | **19** | **20** | **21** | **22** | **23** | **24** |
| f | .11 | .03 | .00 | .03 | .10 | .21 | .34 | .47 | .58 | .68 | .76 | .82 |

Peak is at 15:00. The engine keeps humidity ratio constant over the day (a valid simplification for loads).

---

## 3. Psychrometrics

### 3.1 Equations (ASHRAE Fundamentals Ch.1)

| Property | Equation |
|---|---|
| Atmospheric pressure | `p = 101.325 (1 − 2.25577e-5 · Z)^5.2559`  kPa, Z in m |
| Saturation pressure (t ≥ 0 °C, T in K, Pa) | `ln p_ws = C8/T + C9 + C10·T + C11·T² + C12·T³ + C13·ln T` |
| | C8 = −5.8002206e3, C9 = 1.3914993, C10 = −4.8640239e-2, C11 = 4.1764768e-5, C12 = −1.4452093e-8, C13 = 6.5459673 |
| Saturation over ice (t < 0) | `ln p_ws = C1/T + C2 + C3·T + C4·T² + C5·T³ + C6·T⁴ + C7·ln T` |
| | C1 = −5.6745359e3, C2 = 6.3925247, C3 = −9.6778430e-3, C4 = 6.2215701e-7, C5 = 2.0747825e-9, C6 = −9.4840240e-13, C7 = 4.1635019 |
| Humidity ratio | `W = 0.621945 · p_w / (p − p_w)`  kg/kg dry air |
| Relative humidity | `φ = p_w / p_ws(t)` |
| Enthalpy | `h = 1.006·t + W·(2501 + 1.86·t)`  kJ/kg dry air |
| Specific volume | `v = 0.287042·(t + 273.15)·(1 + 1.607858·W) / p`  m³/kg dry air |
| Density (moist) | `ρ = (1 + W) / v` |
| W from wet bulb t* | `W = [(2501 − 2.326 t*)·W_s* − 1.006 (t − t*)] / (2501 + 1.86 t − 4.186 t*)` |
| Dew point | solve `p_ws(t_d) = p_w` (bisection in engine) |
| Wet bulb | solve `W(t, t*) = W` for t* (bisection) |

### 3.2 Processes

| Process | Rule | Load equation |
|---|---|---|
| Sensible heating/cooling | W constant | `Q = ṁ·c_p·Δt`,  c_p,moist ≈ 1.006 + 1.86 W ≈ 1.02–1.03 kJ/kgK |
| Cooling & dehumidifying | line from entering state to ADP | `Q_total = ṁ (h1 − h2)`; condensate `ṁ (W1 − W2)` |
| Adiabatic mixing | straight line, mass-weighted | `W_m = Σṁ_i W_i / Σṁ_i`, `h_m = Σṁ_i h_i / Σṁ_i` |
| Evaporative (adiabatic) cooling | constant wet bulb | effectiveness `ε = (t1 − t2)/(t1 − t_wb)` = 0.70–0.95 |
| Steam humidification | ~ constant DB | `ṁ_steam = ṁ_a (W2 − W1)` |

### 3.3 Practical airflow constants (sea level, ρ = 1.2 kg/m³)

| Load | Formula (V in L/s) | Formula (V in m³/h) |
|---|---|---|
| Sensible | `Q_s [W] = 1.23 · V · Δt` | `Q_s = 0.34 · V · Δt` |
| Latent | `Q_l [W] = 3010 · V · ΔW` (ΔW kg/kg) | `Q_l = 0.83 · V · ΔW(g/kg)` |
| Total | `Q_t [W] = 1.20 · V · Δh` (kJ/kg) | `Q_t = 0.33 · V · Δh` |

At altitude multiply by `p/101.325`.

### 3.4 Reference states (check your code)

| State | DB | WB | DP | RH | W g/kg | h kJ/kg | v m³/kg |
|---|---|---|---|---|---|---|---|
| Outdoor 36 / 28 WB | 36.00 | 28.00 | 25.49 | 54.9 % | 20.72 | 89.41 | 0.9060 |
| Room 24 °C / 50 % | 24.00 | 17.06 | 12.95 | 50.0 % | 9.31 | 47.84 | 0.8554 |

(p = 101.205 kPa at 10 m.)

---

## 4. Heat transfer through the envelope

### 4.1 U-value

```
R_total = R_si + Σ (d_i / k_i) + Σ R_airgap + R_so          U = 1 / R_total   [W/m²K]
```

| Surface resistance (m²K/W) | Value |
|---|---|
| R_si wall (horizontal heat flow) | 0.13 |
| R_si roof, heat up / heat down | 0.10 / 0.17 |
| R_so (outside, ~4 m/s wind) | 0.04 (h_o ≈ 25); summer design h_o = 17 W/m²K → 0.059 |
| Unventilated air gap 25–100 mm | 0.17–0.18 |

### 4.2 Thermal conductivity of materials (typical, W/mK)

| Material | k | Density kg/m³ |
|---|---|---|
| Clay brick (solid) | 0.72–0.81 | 1800 |
| Hollow clay brick | 0.40–0.55 | 1200 |
| Concrete (dense) | 1.40–1.75 | 2300 |
| AAC block | 0.12–0.16 | 500 |
| Cement plaster | 0.72 | 1860 |
| Gypsum board | 0.16–0.25 | 900 |
| Mineral wool | 0.035–0.040 | 30–100 |
| EPS / XPS | 0.033–0.038 / 0.029–0.034 | 15–35 |
| PU foam | 0.022–0.028 | 35 |
| Glass | 1.0 | 2500 |
| Steel / Aluminium | 50 / 160–230 | 7800 / 2700 |

**Example — 220 mm brick + 2×15 mm plaster:**
R = 0.13 + 0.015/0.72 + 0.22/0.77 + 0.015/0.72 + 0.04 = 0.13 + 0.021 + 0.286 + 0.021 + 0.04 = 0.498 → **U = 2.01 W/m²K**.
Add 50 mm EPS (k = 0.035): R += 1.43 → R = 1.93 → **U = 0.52 W/m²K**.

### 4.3 Typical U, SHGC (glazing whole-window values)

| Glazing | U W/m²K | SHGC | VT |
|---|---|---|---|
| Single clear 6 mm | 5.8 | 0.82 | 0.88 |
| Single tinted 6 mm | 5.8 | 0.60 | 0.50 |
| Single reflective | 5.4 | 0.30–0.45 | 0.20–0.35 |
| Double clear 6-12-6 | 2.8 | 0.70 | 0.78 |
| Double low-e (solar control) | 1.6–1.8 | 0.25–0.40 | 0.40–0.70 |
| Triple low-e argon | 0.6–1.0 | 0.30–0.50 | 0.50–0.70 |
| Frame: aluminium (no break) / thermal break / uPVC | +0.3…+1.0 over centre-glass | | |

Internal shading attenuation (IAC, multiplier on SHGC): light venetian blind 0.6–0.75,
light roller blind 0.4–0.6, drapes 0.5–0.8. External fins/overhang: compute shading geometry.

### 4.4 Dynamic conduction — decrement factor & time lag

Heavy walls store heat: the peak heat flux inside is **reduced (decrement f)** and **delayed (lag φ)**.
The engine uses the CIBSE / admittance formulation with sol-air temperature:

```
q(t) = U·A · [ (T̄_sa − T_i) + f · (T_sa(t − φ) − T̄_sa) ]
```

| Construction (approximate) | f | φ (h) |
|---|---|---|
| Lightweight metal sandwich panel | 0.95–1.0 | 0.5–1 |
| 110 mm brick, plastered | 0.75 | 3 |
| 220 mm brick, plastered | 0.40–0.50 | 6–8 |
| 200 mm hollow concrete block | 0.55 | 5 |
| 150 mm concrete roof + 50 mm insulation | 0.30–0.40 | 5–7 |
| 300 mm concrete | 0.15–0.25 | 9–10 |

For exact values compute from material layers (CIBSE Guide A, ISO 13786) or use ASHRAE CTS tables.

---

## 5. Solar radiation

### 5.1 Solar position (solar time)

```
δ  = 23.45 · sin(360·(284 + n)/365)                    declination, n = day of year
H  = 15 · (AST − 12)                                     hour angle, deg
sin β = cos L·cos δ·cos H + sin L·sin δ                  solar altitude β
sin φ = sin H·cos δ / cos β ;  cos φ = (cos H·cos δ·sin L − sin δ·cos L)/cos β   azimuth from south (+W)
AST = LST + ET/60 + (Long − 15·TZ)/15     (ET = equation of time, min; ±16 min)
```

### 5.2 ASHRAE clear-sky model (Fundamentals 2009+, Ch.14)

```
E0 = 1367 · [1 + 0.033 cos(360 (n − 3)/365)]              extraterrestrial W/m²
m  = 1 / [ sin β + 0.50572 (6.07995 + β)^−1.6364 ]         air mass (β in deg)
ab = 1.454 − 0.406 τb − 0.268 τd + 0.021 τb τd
ad = 0.507 + 0.205 τb − 0.080 τd − 0.190 τb τd
Eb = E0 · exp(−τb · m^ab)     normal beam
Ed = E0 · exp(−τd · m^ad)     horizontal diffuse
```
τb, τd are tabulated monthly per station in ASHRAE climatic data. Typical: τb 0.33–0.60, τd 1.8–2.6
(humid tropics in summer: τb ≈ 0.45–0.60, τd ≈ 2.0–2.3).

### 5.3 Irradiance on a surface (tilt Σ, surface azimuth ψ)

```
γ = φ − ψ
cos θ = cos β · cos γ · sin Σ + sin β · cos Σ
E_beam  = Eb · max(0, cos θ)
Y       = max(0.45, 0.55 + 0.437 cos θ + 0.313 cos²θ)          (vertical surfaces)
E_diff  = Ed · (Y sin Σ + cos Σ)                                  Σ ≤ 90°
E_refl  = (Eb sin β + Ed) · ρ_g · (1 − cos Σ)/2                   ρ_g ≈ 0.2
E_t = E_beam + E_diff + E_refl
```

### 5.4 Sol-air temperature

```
T_sa = T_o + α·E_t / h_o − ε·ΔR / h_o
```
α/h_o: light colour 0.026, dark 0.052 m²K/W (h_o = 17 W/m²K). εΔR/h_o = 3.9 K horizontal, 0 vertical.
Absorptance α: white paint 0.25–0.4, light grey 0.5–0.6, brick 0.65–0.75, dark/black 0.9.

---

## 6. Cooling load calculation

### 6.1 Heat gain ≠ cooling load

*Heat gain* is energy entering the space at an instant. The **radiant** part is first absorbed by
walls/floor/furniture and released later by convection → *cooling load* is smaller and later than
heat gain. The convective part becomes load immediately. Latent gain is always instant.

The engine uses a Radiant-Time-Series-type procedure:

```
Q_load(h) = (1 − F_rad) · q(h) + F_rad · Σ_{j=0..23} r_j · q(h − j)
r_j = (1 − a) · a^j  normalised   (a = 0.3 light, 0.5 medium, 0.7 heavy room)
```
> Note: ASHRAE tabulates RTS factors per construction (Fundamentals Ch.18, Tables 19–20). The
> geometric series above is a simplified, smooth approximation; replace `radiant_time_series()`
> with the tabulated 24 values for certified calculations.

### 6.2 Radiant fractions (ASHRAE Ch.18, Table 14, rounded)

| Gain | Radiant | Convective |
|---|---|---|
| People (sensible) | 0.60 | 0.40 |
| Lighting (LED/fluor., recessed, unvented) | 0.57 (0.48–0.67) | 0.43 |
| Equipment (fan-cooled PCs) | 0.10–0.20 | 0.80–0.90 |
| Equipment (no fan) | 0.50 | 0.50 |
| Wall / roof conduction | 0.46 / 0.60 | 0.54 / 0.40 |
| Window conduction | 0.46 | 0.54 |
| Solar through glass (transmitted) | 1.00 (engine uses 0.90 incl. absorbed part) | — |
| Infiltration, ventilation, latent | 0 | 1.00 |

### 6.3 Load components

| Component | Equation |
|---|---|
| Wall / roof conduction | §4.4 with T_sa |
| Glass conduction | `Q = U·A·(T_o(h) − T_i)` |
| Glass solar | `Q = A · SHGC · IAC · E_t(h)` (shaded area: diffuse only) |
| Partition to unconditioned space | `Q = U·A·(T_adj − T_i)`, T_adj ≈ T_o − 3 K |
| People | `Q_s = N · q_s · CLF/schedule`, `Q_l = N · q_l` |
| Lighting | `Q = W_lpd · A · F_use · F_allowance` (F_al = 1.0 LED, 1.2 fluor. magnetic ballast) |
| Equipment | `Q = W_epd · A · load factor` |
| Infiltration sensible | `Q_s = 1.23 · V(L/s) · (T_o − T_i)` |
| Infiltration latent | `Q_l = 3010 · V(L/s) · (W_o − W_i)` |
| Ventilation (OA) | is a **coil** load, not a room load, unless OA is delivered untreated to the room |

### 6.4 Heat gain from people (ASHRAE Fundamentals Ch.18 Table 1, adjusted, W/person)

| Activity | Typical place | Total | Sensible | Latent |
|---|---|---|---|---|
| Seated at theatre | theatre, cinema | 95 | 65 | 30 |
| Seated, very light work | office, hotel, apartment | 115 | 70 | 45 |
| Moderately active office work | office, hotel | 130 | 75 | 55 |
| Standing, light work; walking | department store, retail | 145 | 75 | 70 |
| Sedentary work | restaurant (incl. food) | 160 | 80 | 80 |
| Light bench work | factory | 235 | 80 | 155 |
| Moderate dancing | dance hall | 265 | 90 | 175 |
| Heavy work | factory, bowling | 440 | 170 | 270 |
| Athletics | gymnasium | 525 | 210 | 315 |

### 6.5 Lighting & equipment densities

| Space | LPD W/m² (ASHRAE 90.1-2019 BAM / typical LED) | Equipment W/m² |
|---|---|---|
| Office | 6.9 / 6–10 | light 5–8, medium 10–15, heavy 15–25 |
| Retail | 9.6–13 / 10–15 | 3–5 |
| School / classroom | 7.4–8.8 | 5–10 |
| Hotel room | 5.4–8 | 3–5 |
| Hospital | 10–12 | 10–20 |
| Server/data room | 5–10 | 500–2000+ (from IT kW) |
| Warehouse | 4–5 | 1–3 |

### 6.6 Infiltration (when building is not pressurised)

| Construction | ACH |
|---|---|
| Tight, pressurised commercial | 0–0.2 |
| Average commercial | 0.3–0.5 |
| Leaky / old, operable windows | 0.5–1.5 |
| Doors to outside: revolving / swing (per use) | add per door flow method |

### 6.7 Safety & diversity

* Room loads: no extra safety (inputs are already design-peak). Optional 5–10 % on fan/coil.
* Diversity on people & lighting for **plant** (block): offices 0.75–0.90.
* Duct heat gain + fan heat: supply fan raises air 0.5–1.0 K per 500 Pa (`ΔT ≈ ΔP / (ρ c_p η)` if motor in stream).

---

## 7. Heating load calculation

Steady state at winter design temperature, no solar or internal credit:

```
Q_heat = [ Σ U_i·A_i + 0.34 · ACH · V_room ] · (T_i − T_o) · (1 + margin)
         0.34 = ρ c_p / 3600  [Wh/m³K]
margin: 10–20 % (pick-up / intermittent heating), + ventilation load at the AHU:
Q_vent = 1.23 · V_OA(L/s) · (T_supply − T_o)  (minus heat recovery: × (1 − ε_HR))
```
Ground floor slab: use ISO 13370 or perimeter method `Q = F_p · P · ΔT` (F_p ≈ 1.4–2.0 W/mK uninsulated).

---

## 8. Ventilation — ASHRAE 62.1

### 8.1 Ventilation Rate Procedure (VRP)

```
V_bz = R_p · P_z + R_a · A_z           breathing-zone outdoor air  [L/s]
V_oz = V_bz / E_z                       zone outdoor air
Single-zone:   V_ot = V_oz
100 % OA:      V_ot = Σ V_oz
Multi-zone recirculating:
   Z_pz = V_oz / V_pz         (V_pz = zone primary supply, use min VAV flow)
   X_s  = V_ou / V_ps         (V_ou = D·ΣR_p P_z + ΣR_a A_z, D = occupant diversity)
   E_v  = 1 + X_s − Z_pz,max  (simplified; or Table 6-3)
   V_ot = V_ou / E_v
```

### 8.2 Table (ASHRAE 62.1-2019/2022 Table 6-1, selected)

| Space | R_p L/s·person | R_a L/s·m² | Default density /100 m² | Combined L/s·person |
|---|---|---|---|---|
| Office space | 2.5 | 0.3 | 5 | 8.5 |
| Conference / meeting | 2.5 | 0.3 | 50 | 3.1 |
| Reception area | 2.5 | 0.3 | 30 | 3.5 |
| Classroom (age 9+) | 3.8 | 0.6 | 35 | 5.4 |
| Lecture classroom | 3.8 | 0.3 | 65 | 4.3 |
| Retail sales | 3.8 | 0.6 | 15 | 7.8 |
| Restaurant dining | 3.8 | 0.9 | 70 | 5.1 |
| Hotel bedroom | 2.5 | 0.3 | 10 | 5.5 |
| Health club / aerobics | 10 | 0.3 | 40 | 10.8 |
| Lobby | 2.5 | 0.3 | 10 | 5.5 |

### 8.3 Zone air-distribution effectiveness E_z

| Configuration | E_z |
|---|---|
| Ceiling supply of cool air | 1.0 |
| Ceiling supply of warm air, floor return | 1.0 |
| Ceiling supply of warm air ≥ 8 K above space, ceiling return | 0.8 |
| Floor supply cool air, ceiling return (UFAD, low velocity displacement) | 1.0–1.2 |

### 8.4 Exhaust rates (62.1 Table 6-4, selected)

| Space | Exhaust |
|---|---|
| Toilet, public | 25 L/s per WC/urinal (or 35 if intermittent) |
| Toilet, private | 12.5 L/s continuous / 25 L/s intermittent |
| Kitchen, commercial | 3.5 L/s·m² (plus hood design) |
| Janitor, trash | 5 L/s·m² |
| Parking garage | 3.7 L/s·m² (or CO control) |

**CO₂ check:** steady state `C_in − C_out = N / V_per_person`; at 2.5 + area component, office indoor CO₂ ≈ 700 ppm above outdoor.
Demand-controlled ventilation (DCV): modulate OA by CO₂ sensor, setpoint ≈ outdoor + 600–700 ppm.

---

## 9. Air-side process: supply air, mixing, coil

### 9.1 Supply airflow

```
V_sa = Q_s,room / (ρ · c_p · (T_room − T_supply))          ρ c_p ≈ 1.2 × 1026 J/m³K
W_sa = W_room − Q_l,room / (ρ · V_sa · 2501000)
RSHR = Q_s / (Q_s + Q_l)
```
Typical supply temperature: 12–14 °C (ΔT 10–12 K) mixing ventilation; 16–18 °C displacement / UFAD;
cold-air distribution 7–10 °C (needs induction diffusers). Airflow sanity check: 4–8 L/s·m² office.

### 9.2 Mixing & coil

```
x_OA = V_OA / V_sa
Mixed:   W_m = x·W_o + (1 − x)·W_r ,  h_m = x·h_o + (1 − x)·h_r
Off-coil state = supply state − fan heat (draw-through) or = supply (blow-through, add fan heat before coil)
ADP: extend line (mixed → off-coil) to the saturation curve
Bypass factor  BF = (t_off − ADP) / (t_on − ADP)       contact factor CF = 1 − BF
Coil total     Q_c = ṁ_a (h_on − h_off)
Coil sensible  Q_cs = ṁ_a c_p,moist (t_on − t_off)
Condensate     ṁ_w = ṁ_a (W_on − W_off)
```
| Coil rows | Typical BF (2.5 m/s face velocity) |
|---|---|
| 4 rows | 0.10–0.20 |
| 6 rows | 0.05–0.10 |
| 8 rows | 0.02–0.05 |

Coil design: face velocity ≤ 2.5 m/s (avoid moisture carry-over), CHW 6/12 °C or 7/13 °C (ΔT 5.5–6 K;
large plants 5/13 °C), air-side ΔP 100–200 Pa wet, water ΔP 30–60 kPa.

When required ADP < ~8 °C or RSHR < 0.65 (humid climates, high OA): use **reheat**,
dedicated outdoor-air system (DOAS), desiccant or heat-pipe/run-around coil.

---

## 10. System selection

| System | Typical use | Capacity range | Pros | Cons |
|---|---|---|---|---|
| Split / multi-split (DX) | small office, residence | 2.5–20 kW | cheap, simple | limited OA, short pipe run |
| VRF / VRV heat pump | offices, hotels | 10–150 kW per outdoor unit | zoning, heat recovery, part-load COP | refrigerant volume (EN 378 check), OA separate |
| Packaged rooftop (DX) | low-rise retail | 10–200 kW | all-in-one, economiser | noise, roof load |
| Chilled water + AHU (CAV) | single large zones: halls, shops | 100 kW–MW | central maintenance, filtration | wasteful for multi-zone |
| Chilled water + VAV | multi-zone offices | MW | fan energy, flexible | humidity at part load, complex controls |
| CHW + FCU + DOAS | hotels, offices, hospitals | MW | decoupled latent/OA, small ducts | condensate drains at each FCU |
| Chilled beams + DOAS | offices, labs | MW | very quiet, low fan energy | dew-point control critical (CHW 14–16 °C) |
| Radiant slab + DOAS | high-end, airports | MW | comfort, low exergy | slow response |
| District cooling / thermal storage | campuses | MW | peak shaving | capex |

---

## 11. Air distribution & diffusers

| Parameter | Recommended |
|---|---|
| Throw: terminal velocity | 0.25 m/s at the wall / 0.75 × room length |
| Occupied zone velocity | ≤ 0.20 m/s (cooling) |
| ADPI (Air Diffusion Performance Index) | ≥ 80 % |
| Diffuser neck velocity | 2–3 m/s offices (NC 30–35) |
| Return grille face velocity | 1.5–2.5 m/s (occupied), up to 4 m/s high level |
| Diffuser ΔP total | 15–30 Pa |
| Ceiling diffuser spacing | ≈ 1.5–2 × ceiling height, square diffuser 3×3 to 4×4 m grid |

Noise criteria: office NC 35–40, meeting NC 30, hotel room NC 25–30, studio NC 15–20.

---

## 12. Duct design

### 12.1 Pressure loss

```
Friction (Darcy–Weisbach):   Δp_f = f · (L/D) · ρ V²/2
Colebrook:   1/√f = −2 log10( ε/(3.7 D) + 2.51/(Re √f) ) ,  Re = V D / ν   (ν_air = 1.53e-5 m²/s)
Fitting:     Δp = C · ρ V²/2          (C from ASHRAE Duct Fitting Database / SMACNA)
Velocity pressure p_v = 0.602 V²  [Pa]  (ρ = 1.204)
Total system:  Δp_total = Σ (friction + fittings) along the critical (index) path + components
```

Roughness ε: galvanised steel 0.09 mm, aluminium 0.05 mm, fibreglass duct-board 0.9 mm,
flexible duct (fully extended) 3 mm (avoid > 1.5 m runs).

### 12.2 Equivalent diameter (rectangular → round, Huebscher)

```
D_e = 1.30 · (a·b)^0.625 / (a + b)^0.25
```
Keep aspect ratio ≤ 4:1 (max 8:1). Flat-oval: `D_e = 1.55 A^0.625 / P^0.25`.

### 12.3 Design velocities & friction (commercial, low pressure)

| Section | Velocity m/s (office, NC 35) | Velocity m/s (industrial) |
|---|---|---|
| Main supply duct | 5–8 | 10–12 |
| Branch duct | 3–5 | 6–9 |
| Run-out to diffuser | 2–3 | 4–5 |
| Return main | 4–6 | 8–10 |
| AHU coil face | 2.0–2.5 | 2.5–3.0 |
| Filter face | 1.5–2.5 | 2.5 |
| Louvres (intake, net free area) | 2.0–2.5 | 3.5 |

| Method | Rule | Use |
|---|---|---|
| Equal friction | constant 0.8–1.2 Pa/m (low pressure) | most comfort systems |
| Static regain | size so static regain offsets the loss in next section | high-velocity / VAV mains |
| Velocity reduction | choose velocities step by step | simple systems |
| T-method | optimisation (life-cycle cost) | large systems |

### 12.4 Common fitting coefficients C (approximate, based on downstream/section velocity)

| Fitting | C |
|---|---|
| 90° elbow, round, smooth radius R/D = 1.5 | 0.15–0.25 |
| 90° elbow, rectangular, mitred with turning vanes | 0.25–0.35 |
| 90° elbow, mitred, no vanes | 1.1–1.3 |
| Branch tee (branch side, 90°) | 0.6–1.5 |
| Tee straight-through | 0.05–0.3 |
| Sudden expansion (A1/A2 = 0.5) | 0.25 (upstream V) |
| Contraction 30° | 0.05 |
| Volume damper fully open | 0.2–0.5 |
| Duct entry from plenum (sharp) / bellmouth | 0.5 / 0.03 |
| Exit to room | 1.0 |

### 12.5 Equal-friction procedure (engine `size_duct_network`)

1. Draw single-line diagram; give each section a name, flow, length, upstream section.
2. Estimate ΣC per section (elbows, tees, dampers, transitions).
3. For each section solve the exact D giving R_target (e.g. 1.0 Pa/m); also check `D ≥ √(4Q/πV_max)`.
4. Round up to the standard size (EN 1506 series: 100, 125, 160, 200, 250, 315, 355, 400, 450, 500, 560, 630, 710, 800, 900, 1000, 1120, 1250 mm).
5. Convert to rectangular for ceiling depth (e.g. 250 mm high) with Huebscher.
6. Compute Δp for each section = R·L + ΣC·p_v, then sum each root→leaf path; the largest = index path.
7. Balance other paths with dampers (or resize to bleed excess pressure).
8. External static pressure (ESP) = index path + terminal (diffuser) losses; fan total = ESP + AHU internal (coil, filters, sections).

### 12.6 AHU internal pressure drops (typical)

| Component | Clean Pa | Design Pa |
|---|---|---|
| Panel filter G4 (MERV 8) | 50 | 150 |
| Bag filter F7 (MERV 13) | 100 | 250 (final 450) |
| HEPA H13 | 250 | 500 |
| Cooling coil, wet, 6 row | 120–180 | |
| Heating coil 2 row | 40–60 | |
| Heat wheel | 100–200 | |
| Silencer | 30–80 | |
| Mixing box, dampers | 30–50 | |

---

## 13. Hydronic (pipe) design

### 13.1 Flow and loss

```
V̇ = Q / (ρ · c_p · ΔT)          c_p,water = 4.19 kJ/kgK,  1 kW at 5.5 K → 0.0434 L/s
Δp = f (L/D) ρV²/2 + Σ K ρV²/2         (same Colebrook, steel ε = 0.045 mm, copper 0.0015 mm)
Water ρ (Thiesen): 1000·[1 − (t + 288.9414)/(508929.2·(t + 68.12963))·(t − 3.9863)²]
Water μ (Vogel):   2.414e-5 · 10^(247.8/(T − 140))   Pa·s, T in K
```

### 13.2 Design criteria

| Parameter | Value |
|---|---|
| Friction rate | 100–400 Pa/m (common design 250 Pa/m); ≤ 300 Pa/m for quiet areas |
| Velocity (≤ DN50) | 0.6–1.2 m/s |
| Velocity (DN65–DN150) | 1.2–2.4 m/s |
| Velocity (> DN150) | 2.0–3.0 m/s (max 3.0 m/s erosion/noise) |
| CHW ΔT | 5.5–6 K (standard), 8 K (low-flow/large plants) |
| Condenser water ΔT | 5–5.6 K (flow ≈ 0.054 L/s per kW cooling) |
| Heating hot water | 80/60 °C classic, 60/40 or 45/35 °C condensing boilers/heat pumps |
| Fittings allowance (first estimate) | equivalent length = 50–100 % of straight length |

### 13.3 Steel pipe Sch.40 internal diameters

| DN | 15 | 20 | 25 | 32 | 40 | 50 | 65 | 80 | 100 | 125 | 150 | 200 | 250 | 300 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ID mm | 15.8 | 20.9 | 26.6 | 35.1 | 40.9 | 52.5 | 62.7 | 77.9 | 102.3 | 128.2 | 154.1 | 202.7 | 254.5 | 303.2 |
| Max L/s @ 250 Pa/m | 0.09 | 0.19 | 0.36 | 0.76 | 1.15 | 2.25 | 3.61 | 6.45 | 13.3 | 24.2 | 39.4 | 81.2 | 148 | 234* |
| Velocity m/s at that flow | 0.44 | 0.55 | 0.65 | 0.79 | 0.88 | 1.04 | 1.17 | 1.35 | 1.62 | 1.88 | 2.11 | 2.52 | 2.91 | 3.24* |

(Computed with the engine: water 9 °C, ε = 0.045 mm. \*DN300 exceeds 3 m/s at 250 Pa/m — limit by velocity, ≈ 216 L/s.)

### 13.4 Minor-loss coefficients K (water)

| Fitting | K |
|---|---|
| Elbow 90° standard / long radius | 0.75 / 0.45 |
| Tee through / branch | 0.4 / 1.0–1.8 |
| Gate valve open / butterfly open | 0.15 / 0.5–1.0 |
| Globe valve open | 6–10 |
| Check valve swing | 2.0 |
| Strainer (clean) | 2–5 (or use manufacturer ΔP) |

### 13.5 Pump head build-up (index circuit)

```
H = Σ pipe friction (supply + return, index run) + fittings
  + chiller evaporator ΔP (40–80 kPa) + coil ΔP (30–60 kPa) + control valve ΔP (authority ≥ 0.5, ≈ coil ΔP)
  + balancing valve (3–10 kPa) + strainer (10–20 kPa) + 10–15 % margin
```
Control-valve authority `N = Δp_valve / (Δp_valve + Δp_circuit)` ≥ 0.25–0.5.
Valve sizing: `Kv = V̇[m³/h] / √(Δp[bar])`.

### 13.6 System arrangements

| Arrangement | Description |
|---|---|
| Constant primary | simple, 3-way valves, high pump energy |
| Primary–secondary (decoupled) | constant primary per chiller, variable secondary (2-way valves, VSD) |
| Variable primary flow (VPF) | VSD primary pumps, min-flow bypass; lowest energy, needs good controls |
| Expansion tank | acceptance volume = system volume × expansion (≈ 0.5 % for 6→35 °C CHW; 3–4 % for 10→80 °C HW) |

---

## 14. Fans & pumps

```
Air power         P_air = V̇ · Δp_t
Electric input    P_e = V̇ · Δp_t / (η_fan · η_motor · η_drive)
SFP (specific fan power) = P_e / V̇  [W/(m³/s)] ;  target ≤ 1500–2000 (supply+extract AHU ≤ 2500)
Pump input        P_e = V̇ · Δp / (η_pump · η_motor)
Affinity laws (same fan/pump, speed n):  V̇ ∝ n,  Δp ∝ n²,  P ∝ n³
Fan heat to air   ΔT = Δp / (ρ c_p · η_overall)   (motor in air stream)
```

| Fan type | η peak | Use |
|---|---|---|
| Forward-curved centrifugal | 0.55–0.70 | small FCU, low pressure |
| Backward-curved / airfoil centrifugal (plug/EC) | 0.70–0.85 | AHU |
| Axial (vane) | 0.70–0.85 | high flow, tunnels, towers |
| Mixed flow | 0.70–0.80 | inline duct |

Motor efficiency IE3: 7.5 kW ≈ 0.91, 30 kW ≈ 0.935, 110 kW ≈ 0.955. VSD loss ≈ 2–4 %.
Select fan near peak efficiency, right of surge/stall line; system curve `Δp = k·V̇²`.

---

## 15. Central plant: chillers, towers, boilers

### 15.1 Sizing

```
Chiller capacity  = coincident block coil load × (1.05–1.15)   [kW]   (1 TR = 3.517 kW)
Number            = N+1 redundancy for critical; 2–3 equal machines (or 1 small + 2 large) for part-load
Heat rejection    = Q_evap · (1 + 1/COP)
CHW flow          = Q / (ρ c_p ΔT)
Condenser flow    = Q_rej / (ρ c_p ΔT_cw)
Tower make-up     = evaporation (≈ 0.00153 × flow × ΔT/°C in fraction) + drift (0.001–0.005 %) + blowdown
   evaporation ≈ Q_rej / h_fg  (≈ 1.5 L/h per kW rejected)
   blowdown    = evaporation / (cycles − 1)   (cycles of concentration 3–6)
```

### 15.2 Typical performance (full load, AHRI conditions)

| Equipment | COP (kW/kW) | kW/TR |
|---|---|---|
| Split DX inverter | 3.2–4.5 (SEER 5–7) | 0.8–1.1 |
| VRF | 3.5–4.5 (IPLV 6–8) | 0.8–1.0 |
| Air-cooled scroll/screw chiller | 2.9–3.3 | 1.1–1.2 |
| Water-cooled screw chiller | 5.0–5.8 | 0.6–0.7 |
| Water-cooled centrifugal (VSD) | 6.1–7.0 | 0.50–0.58 |
| Absorption double-effect | 1.2–1.4 (thermal) | — |
| Condensing boiler | 0.92–0.98 (GCV) | — |
| Air-source heat pump heating | 2.5–4.0 | — |

Cooling tower: approach 3–4 K (CW leaving − ambient WB), range 5–5.6 K, fan ≈ 0.01–0.02 kW per kW rejected.
Overall plant efficiency (chiller + pumps + towers): good 0.75–0.85 kW/TR, excellent ≤ 0.65 kW/TR.

---

## 16. Energy, efficiency & controls

| Measure | Typical saving |
|---|---|
| VSD on fans & pumps (VAV, VPF) | 30–60 % of fan/pump energy |
| Supply air temperature reset (12→16 °C at low load) | 5–15 % |
| CHW temperature reset / higher CHW (7→9 °C) | 1.5–3 % chiller energy per K |
| Condenser water reset (lower CW when WB low) | 1–3 % per K |
| Air-side economiser (dry climates) | 10–40 % |
| DCV with CO₂ | 10–30 % ventilation energy |
| Energy recovery wheel (ε 0.70–0.80 sensible+latent) | large in humid climates (OA load ↓ 60–75 %) |
| Thermal storage (ice/CHW) | peak demand shift |

**Core control loops (PID, see the Control field for tuning):**

| Loop | Sensor | Actuator | Typical setpoint |
|---|---|---|---|
| Supply air temp | SAT sensor | CHW valve | 13 °C (reset 12–16) |
| Duct static pressure | Static sensor 2/3 down main | Supply fan VSD | 250–375 Pa (trim & respond reset) |
| Zone temp (VAV) | Zone sensor | VAV damper (+ reheat) | 24 ±1 °C |
| Min OA | Airflow station / CO₂ | OA damper | V_ot / DCV |
| Building pressure | ΔP indoor–outdoor | Return/relief fan | +5…+12 Pa |
| CHW ΔP | Remote ΔP sensor | Secondary pump VSD | 60–100 kPa |
| CHW supply temp | Chiller leaving sensor | Chiller (internal) | 6–7 °C |
| CW temp | Tower basin sensor | Tower fan VSD | 29 °C or WB + approach |
| Min chiller flow | Flow meter | Bypass valve | per chiller spec |

---

## 17. Worked example (full output)

**Project:** open-plan office zone, 150 m² × 3.0 m, top floor, west & south façades, latitude 21 °N,
July 21 (n = 202). Design 36 °C DB / 28 °C MCWB, DR = 7 K. Indoor 24 °C / 50 %.

| Element | Data |
|---|---|
| West wall | 18 m² net, U 0.60, α 0.6, f 0.45, φ 7 h |
| South wall | 27 m² net, U 0.60, α 0.6, f 0.45, φ 7 h |
| Roof | 150 m², U 0.45, α 0.7, f 0.35, φ 6 h |
| West glass | 12 m², U 2.8, SHGC 0.40, blinds IAC 0.8 |
| South glass | 18 m², U 2.8, SHGC 0.40, blinds IAC 0.8 |
| People | 15 × (75 S + 55 L) W, 08–18 h |
| Lighting / Equipment | 9 W/m² / 12 W/m² (15 % night base) |
| Infiltration | 0.3 ACH; room = medium weight |

Run: `python hvac_engine.py`. Result (selected):

```
[1] Outdoor : DB 36.00 | WB 28.00 | DP 25.49 | RH 54.9 % | W 20.72 g/kg | h 89.41 kJ/kg
    Room    : DB 24.00 | WB 17.06 | DP 12.95 | RH 50.0 % | W  9.31 g/kg | h 47.84 kJ/kg

[2] hour  T_out  walls  roof  g_cnd  solar  peopl  light  equip  infil | SENS   LAT  TOTAL
      15   36.0    276  1290    989   2893   1120   1344   1798    534 | 10244  2087  12331
      16   35.8    299  1492    985   3174   1122   1347   1799    525 | 10743  2087  12830
      17   35.3    323  1662    956   3103   1124   1348   1799    503 | 10819  2087  12907
    PEAK at hour 17: sensible 10.82 kW, latent 2.09 kW, total 12.91 kW, RSHR 0.84, 86 W/m2

[3] Vbz = 2.5*15 + 0.3*150 = 82.5 L/s

[4] Supply airflow 812 L/s (5.4 L/s.m2, 6.5 ACH), OA fraction 10.2 %
    Supply 13.0 C / 8.44 g/kg ; Mixed 25.24 C / 10.47 g/kg ; Off-coil 12.2 C (fan heat 0.8 K)
    ADP 11.20 C, bypass factor 0.071
    Coil load: total 17.46 kW, sensible 12.63 kW, condensate 6.90 kg/h, SHR 0.72

[5] section   flow L/s   D mm   rect @250mm   V m/s   R Pa/m   dP Pa
    AHU-A         812    450     750x250     5.11    0.62    23.7
    A-B           406    315     350x250     5.21    0.99    10.9
    B-D1/B-D2     203    250     250x250     4.14    0.86    18.9 / 24.5
    A-C           406    315     350x250     5.21    0.99    19.7
    C-D3/C-D4     203    250     250x250     4.14    0.86    18.0 / 23.7
    Critical path C-D4: 67.1 Pa ; fan total ≈ 412 Pa → 536 W input

[6] Heating (22 / 8 °C, 0.5 ACH, +15 %): 4106 W

[7] Plant for 40 zones: chiller 768 kW (218 TR) as 2 × 384 kW, compressor 140 kW,
    tower rejection 908 kW, CHW 33.3 L/s → DN150 (1.79 m/s, 183 Pa/m), CW 43.6 L/s,
    CHW pump at 250 kPa ≈ 12.0 kW
```

**Hand check of the coil load:**
Room 12.91 kW + OA `1.2 × 82.5 × (89.41 − 47.84)/1000 = 4.12 kW` + fan heat `1.23 × 812 × 0.8/1000 = 0.80 kW`
= 17.8 kW ≈ 17.46 kW ✔ (difference = density and mass-weighting).

**Interpretation:**
* Peak at 17:00 because of the west glass (solar 3.1 kW = 29 % of sensible) and the roof lag.
  → External shading on the west façade is the biggest lever.
* 86 W/m² is typical for a top-floor office in a hot-humid climate (range 80–150 W/m²).
* ADP 11.2 °C with BF 0.07. The coil rating in §C1.4 selects **8 rows × 10 FPI** (3/8" tube, 900 × 381 mm face, 2.37 m/s) with 6/12 °C CHW.

---

## 18. Engine API reference

```python
from hvac_engine import *

# Psychrometrics
p  = pressure_at_altitude(500)                 # kPa
s  = AirState.from_rh(24, 0.5, p)              # or AirState.from_wb(36, 28, p) / AirState(t, w, p)
s.h, s.w, s.rh, s.twb, s.tdp, s.v, s.rho
m  = mix(stateA, massA, stateB, massB)

# Loads
site = Site(latitude=10.8, longitude=106.7, tz=7, tau_b=0.5, tau_d=2.1)
zone = Zone("Z1", floor_area=100, height=3, people=10)
zone.walls.append(Opaque("E", 20, 0.6, tilt=90, azimuth=-90, absorptance=0.6, decrement=0.5, lag_h=6))
zone.windows.append(Window("E-g", 8, 2.8, 0.4, azimuth=-90))
r = cooling_load(zone, site, day=172, t_db_max=35, daily_range=8, t_wb_mc=27)
r["sensible"][h], r["latent"][h], r["peak_hour"], r["components"]["solar"]
heating_load([(U, A), ...], t_in=22, t_out=5, volume=300, ach=0.5)

# Ventilation
Vbz, Voz = outdoor_air("classroom_9plus", area=70, people=30, Ez=1.0)
Vot, Ev  = system_oa_multizone([(Voz1, Vpz1), (Voz2, Vpz2)])

# Air process
V   = supply_airflow(q_sens, t_room=24, t_supply=13)      # m3/s
w_s = supply_humidity(room_state, q_lat, V)
adp, bf = coil_adp_bf(mixed_state, off_coil_state)
Qt, Qs, condensate = coil_load(mixed_state, off_coil_state, V)

# Ducts
D   = size_round_for_friction(0.5, 1.0)                    # m
De  = equiv_diameter_rect(0.6, 0.3)
rows, paths = size_duct_network([DuctSection("A", 1.0, 10, 0.5), DuctSection("B", 0.5, 8, 1.2, "A")])

# Pipes & plant
q   = water_flow(500e3, 5.5)                               # m3/s
size_pipe(q, max_pa_m=300)                                 # {'DN':..,'V':..,'R':..}
plant_sizing(coil_total_w, cop=6.0, n_chillers=3)
fan_power(V, dp_pa), pump_power(q, head_pa), fan_laws(Q1, dp1, P1, n1, n2)

# Equipment & BOM (procurement)
motor_size(5.6)                                            # -> 7.5 kW IEC
control_valve_kv(2.5, 40)                                  # -> (Kv_calc, Kvs_std, dp_actual_kPa)
select_diffusers(0.812, 150)                               # -> (count, m3/s each, neck mm)
duct_sheet_metal(rows)                                     # -> [(name, D, m2, t_mm, kg)]
expansion_tank(6000, static_bar=3.0)                       # -> (Ve L, Vtank L)
items = zone_bom(...) + plant_bom(coil_total_w, n_zones=40)
print_bom(items); bom_to_csv(items, "bom.csv")             # Excel-ready, add unit prices

# Equipment engineering (Part C)
g  = CoilGeometry(face_w=0.9, face_h=0.381, rows=8, fpi=10, circuits=11)
r  = rate_coil(g, entering_state, 0.812, tw_in=6.0, q_water=0.00069)   # Q, SHR, air_out, dp...
g, r = design_coil(entering_state, 0.812, Q_req=17460, t_leave_max=12.5)
lmtd(dt1, dt2); eff_counterflow(ntu, cr); eff_crossflow_unmixed(ntu, cr); nu_tube(Re, Pr)
fin_efficiency_schmidt(h, k, t, r_o, Pt, Pl)
insulation_condensation(6, 32, 0.80, d_out_pipe=0.0603)   # -> (thickness m, dew point)
npsh_available(101.3, 32, z_static_m=1.5, h_friction_m=0.8)
vapour_compression(h1, h2s, h3, Q_evap_kw, eta_is=0.72, v1=0.0709)
db_sum([75, 70, 70]); room_spl(Lw, r, R); vibration_isolation(1450, 25)
duct_leakage_en12237(area_m2, p_pa, 'B')
```

**Extending the engine:**
* Replace `radiant_time_series()` with ASHRAE tabulated RTS for your construction.
* Add `Opaque` interior partitions with `tilt=90` and set an adjacent-temperature model.
* Loop `cooling_load` over 12 months (change `day`, `tau_b`, `tau_d`, design temps) to find the true annual peak (south façades in the northern tropics often peak in winter months!).
* Sum `r["sensible"]` hour-by-hour over zones → block load for AHU/plant.

---

## 19. Rules of thumb (for checking only)

| Building | Cooling W/m² (hot-humid) | Supply air L/s·m² | CHW L/s per kW |
|---|---|---|---|
| Office | 100–150 (perimeter) / 70–100 (core) | 4–8 | 0.043 (ΔT 5.5) |
| Retail / mall | 150–250 | 6–10 | |
| Restaurant | 200–350 | 8–12 | |
| Hotel room | 70–110 | 3–5 | |
| Hospital | 120–200 | 6–12 | |
| Auditorium | 200–350 (≈ 0.15 kW/seat) | 8–15 | |
| Data centre | 500–3000 | — | |

1 TR = 3.517 kW ≈ 12 000 BTU/h · 1 kW ≈ 3412 BTU/h · 1 L/s = 2.119 CFM = 3.6 m³/h · 1 Pa = 0.004 in.wg.
Typical office: ≈ 10–12 m²/kW (≈ 30–40 m²/TR).

---

## 20. Design checklist

- [ ] Zoning by orientation, use, schedule; interior vs perimeter separated
- [ ] Climate data source & percentile documented; dew-point design checked for humid climates
- [ ] U-values meet energy code (e.g. QCVN 09, ASHRAE 90.1 prescriptive tables)
- [ ] Hourly loads for several months (not only July) — find each zone's true peak
- [ ] OA per 62.1 (or local), exhaust balanced, building positively pressurised (+5–12 Pa)
- [ ] Psychrometric process plotted: ADP achievable (≥ 8 °C with 6 °C CHW), RH in room ≤ 60 %
- [ ] Supply ΔT and airflow consistent with diffuser selection (ADPI ≥ 80 %, NC target)
- [ ] Ducts: velocities & noise; index path ESP; dampers on branches; fire dampers at compartment walls
- [ ] Pipes: velocity/friction; index circuit pump head; control valve authority; air vents/drains; expansion
- [ ] Plant: block load, part-load strategy, N+1, COP/IPLV vs code minimum
- [ ] Condensate drains with traps (coil negative pressure: trap height ≥ ESP in mm + 25 mm)
- [ ] Controls sequence written; BMS points list; commissioning & TAB (testing, adjusting, balancing) plan
- [ ] Space for maintenance (coil pull, filter access), acoustic & vibration isolation

---

## 21. Equipment selection & bill of materials (procurement)

Design gives **capacities**; procurement needs **specifications** a supplier can quote against.
For every item send: duty (capacity, flow, pressure), design conditions, electrical supply, limits
(noise, efficiency, dimensions), standards and accessories. Compare offers on the **same duty point**.

### 21.1 From design result to purchase spec

| Item | Selection rule | Spec to send to supplier | Engine |
|---|---|---|---|
| AHU / ducted FCU | airflow, coil kW (total + sensible), ESP | L/s, ESP Pa, coil EAT/LAT (DB/WB), CHW in/out °C and L/s, rows, face velocity ≤ 2.5 m/s, filter classes, motor type (EC/VSD), casing (double-skin 25–50 mm), drain pan SS | `zone_bom` |
| Ceiling diffuser | one per 9–20 m² (≈ 4×4 m), ≤ 150 L/s each, neck velocity ≤ 3 m/s | neck Ø, L/s, throw at 0.25 m/s, NC ≤ 35, plenum box + damper | `select_diffusers` |
| Return grille | face velocity ≤ 2.0 m/s | size, L/s, filter, finish | `select_return_grilles` |
| Duct sheet metal | area = πDL × 1.15; thickness by largest side (below) | Z275 galvanised, thickness, kg, sealing class | `duct_sheet_metal` |
| Control valve | Kv = Q[m³/h]/√Δp[bar], Δp_valve ≈ Δp_coil (authority ≈ 0.5) | DN, Kvs, equal-% characteristic, actuator 24 VAC 0–10 V, close-off pressure; or PICV with flow range | `control_valve_kv` |
| Pump | flow, head; motor = next IEC size ≥ shaft × 1.15 | L/s, kPa, NPSH available, efficiency ≥ 75 %, IE3 motor, VSD, N+1 | `plant_bom`, `motor_size` |
| Chiller | block load × 1.1, split into 2–3 units | kW each, CHW & CW temperatures/flows, COP and IPLV, refrigerant, max water Δp, voltage, sound dB(A) | `plant_bom` |
| Cooling tower | Q_rej = Q(1 + 1/COP) | kW, L/s, in/out °C, design WB, approach, drift, fan VSD, materials | `plant_bom` |
| Expansion tank | V_t = V_e / (1 − (P_i+1)/(P_f+1)) | litres, PN, pre-charge bar | `expansion_tank` |
| Insulation | prevent condensation (surface > ambient dew point) | material, thickness, density, vapour barrier | `chw_insulation` |

**Galvanised duct thickness (low pressure, ≤ 500 Pa):**

| Largest side / Ø (mm) | ≤ 300 | 301–750 | 751–1500 | 1501–2100 | > 2100 |
|---|---|---|---|---|---|
| Sheet (mm) | 0.5 | 0.6 | 0.8 | 1.0 | 1.2 |

Mass = area × t × 7850 kg/m³ (0.6 mm ≈ 4.7 kg/m²). Supports every ≤ 2.4 m (M8 rod for ducts ≤ 1000 mm, M10 above).

**CHW pipe insulation (closed-cell elastomeric, high-ε surface, hot-humid ambient 32 °C / 80 % RH), derived in §C6:**

| Pipe | ≤ DN50 | DN65–150 | DN200–300 | > DN300 |
|---|---|---|---|---|
| Thickness (mm) | 25 | 32 | 40 | 50 |

With bright aluminium/SS jackets, add one size (see §C6).

Duct insulation: glass wool 25 mm, 24–32 kg/m³ + FSK (supply ducts in ceiling); 50 mm outdoors with cladding;
or elastomeric 19–25 mm.

**Standard Kvs series:** 0.25, 0.4, 0.63, 1.0, 1.6, 2.5, 4.0, 6.3, 10, 16, 25, 40, 63, 100, 160, 250, 400.
**IEC motor sizes (kW):** 0.37, 0.55, 0.75, 1.1, 1.5, 2.2, 3, 4, 5.5, 7.5, 11, 15, 18.5, 22, 30, 37, 45, 55, 75, 90, 110, 132, 160, 200.

### 21.2 Worked-example equipment schedule (from `python hvac_engine.py`)

**Zone (one of 40), 150 m² office:**

| Qty | Item | Specification |
|---|---|---|
| 1 set | AHU / ducted FCU | 812 L/s (2925 m³/h); ESP 92 Pa; fan TSP 412 Pa; EC motor 0.75 kW; 6-row CHW coil, face 0.32 m²; total 17.5 kW / sensible 12.6 kW; EAT 25.2/18.4 WB, LAT 12.2/11.8 WB; CHW 6/12 °C 0.69 L/s; G4 + F7 filters |
| 10 pcs | Square ceiling diffuser 600×600, 4-way | neck Ø200, 81 L/s each, NC ≤ 35, plenum box + damper |
| 2 pcs | Return grille 600×600 | face velocity ≤ 2.0 m/s |
| 1 set | OA louvre + motorised damper | 82 L/s, OA duct Ø200 |
| 70.9 kg | Galvanised sheet 0.5 mm | Ø250 branches (18.1 m²) |
| 147 kg | Galvanised sheet 0.6 mm | Ø450 and Ø315 mains (31.2 m²) |
| 47.1 m² | Duct insulation | glass wool 25 mm + FSK |
| 15 m | Flexible insulated duct Ø200 | ≤ 1.5 m per diffuser |
| 6 pcs | Volume control dampers | one per branch |
| 19 sets | Duct supports | rod + angle, ≤ 2.4 m spacing |
| 30 m | CHW pipe DN32 | 0.72 m/s, 211 Pa/m (supply + return) |
| 30 m | Elastomeric insulation 25 mm | for DN32 |
| 15 m | Condensate uPVC DN25–32 + P-trap | 6.9 kg/h, trap seal ≥ 53 mm |
| 1 set | 2-way control valve DN32 + actuator | Kvs 4.0, Δp 39 kPa, authority 0.49 (or PICV set to 0.69 L/s) |
| 2 / 1 / 1 / 2 | Isolating valves / Y-strainer / balancing valve / flex connectors | DN32 PN16 |
| 2 / 2 / 1 / 1 | Pressure gauges / thermometers / air vent / drain valve | |
| 1 each | Zone T+RH sensor, CO₂ sensor, SAT sensor, filter ΔP switch, fan status, drain-pan float switch, DDC controller (≈ 8 AI, 6 AO, 6 DI, 4 DO) | |

**Central plant (40 zones, block load 768 kW):**

| Qty | Item | Specification |
|---|---|---|
| 2 | Water-cooled screw chiller | 384 kW (109 TR) each; CHW 12→6 °C 16.7 L/s; CW 32→37 °C 21.8 L/s; COP ≥ 5.5, IPLV ≥ 6.9; R134a / R513A / R1234ze |
| 2 | Cooling tower, induced-draft counter-flow FRP | 454 kW each; 21.8 L/s; 37→32 °C at 28 °C WB; VSD fan |
| 2+1 | CHW pumps | 16.7 L/s (60 m³/h) @ 250 kPa (25 m); 7.5 kW IE3 + VSD |
| 2+1 | CW pumps | 21.8 L/s (78 m³/h) @ 200 kPa (20 m); 7.5 kW IE3 + VSD |
| — | CHW / CW headers | DN150 (1.79 m/s) / DN150 (2.34 m/s); length from layout |
| 1 | Diaphragm expansion tank | system ≈ 6150 L (8 L/kW), expansion 37 L, tank ≥ 111 L PN10, pre-charge 3.3 bar |
| 1 lot | Air/dirt separator, chemical dosing, make-up tank (≈ 1360 L/h evaporation), BTU meter, plant manager, remote ΔP sensor, plant valves | |

`bom_to_csv()` writes `bom_example.csv` with `unit_price` and `amount = qty × price` columns. Fill in
supplier prices and the totals calculate in Excel or LibreOffice. Prices are not included on purpose:
they change with market, region and time, so always use current quotations.

### 21.3 Market categories & example manufacturers

These are well-known examples only, not recommendations. Always compare at least 3 offers on the same duty point.

| Category | Example brands (international / widely available in Asia) |
|---|---|
| Chillers | Daikin, Carrier, Trane, York (Johnson Controls), Mitsubishi Heavy, LG, Gree, Midea |
| AHU / FCU | Daikin, Carrier, Trane, Systemair, TROX, Kruger, local OEM assemblers |
| VRF / split | Daikin, Mitsubishi Electric, LG, Samsung, Panasonic, Toshiba, Hitachi |
| Cooling towers | BAC, Evapco, Marley (SPX), Liang Chi, Kuken, Tashin |
| Pumps | Grundfos, Wilo, KSB, Ebara, Xylem (Lowara, Bell & Gossett) |
| Valves, PICV, actuators | Belimo, Danfoss, Siemens, Honeywell, IMI Hydronic (TA), Oventrop |
| Diffusers, grilles, dampers | TROX, Price, Titus, Krueger, local fabricators |
| Fans | Ziehl-Abegg, ebm-papst, Systemair, Kruger, Nicotra-Gebhardt |
| VSD | ABB, Danfoss, Schneider, Siemens, Yaskawa, Mitsubishi |
| BMS / DDC | Siemens, Honeywell, Schneider, Johnson Controls, Distech, Delta Controls |
| Sensors | Siemens, Belimo, Vaisala, E+E, Dwyer, Sensirion (OEM CO₂/RH) |
| Insulation | Armacell (Armaflex), K-Flex, Kingspan, Owens Corning, Rockwool |
| Pipe & fittings | steel ASTM A53/A106 Sch.40; grooved: Victaulic, Shurjoint; PP-R for small CHW |

### 21.4 Procurement checklist

- [ ] Data sheet per item, with duty point, accessories and certification (AHRI, Eurovent, CE/UL)
- [ ] Supplier selection-software printout showing performance at **your** conditions (not catalogue nominal)
- [ ] Part-load data (chillers: IPLV/NPLV; pumps and fans: curves plotted with the system curve)
- [ ] Electrical: voltage, FLA, starting method, power factor, IP rating
- [ ] Sound data (dB(A) at 1 m, or octave-band Lw) for AHUs, chillers and towers
- [ ] Dimensions, weight, access and lifting points; delivery route through doors and shafts
- [ ] Spares, warranty (≥ 1–2 years, compressor 5 years), local service, lead time
- [ ] Factory acceptance test (FAT) for chillers and AHUs on large jobs

---

## 22. Installation tools & test (TAB) instruments

### 22.1 Fabrication & installation

| Work | Tools / machines |
|---|---|
| Duct fabrication | CNC plasma cutter, Pittsburgh lock former, TDF/TDC flange roll former, folding brake, guillotine shear, spiral duct machine (round), beading machine, hand seamer, sealant gun |
| Duct installation | rotary hammer drill, anchors, threaded-rod cutter, scissor lift, rivet gun, flange clamps, duct sealant + foil tape |
| Pipe (steel) | band saw / pipe cutter, threading machine (≤ DN50), roll-groover (grooved joints), SMAW/TIG welding set, bevelling machine, chain blocks |
| Pipe (copper / refrigerant) | tube cutter, flaring and swaging tools, bender, brazing torch (oxy-acetylene) with 5–15 % Ag rods, nitrogen cylinder + regulator (purge while brazing) |
| Refrigerant (DX / VRF) | vacuum pump (≥ 6 CFM, to ≤ 500 micron), micron gauge, digital manifold gauges, refrigerant scale, recovery unit, electronic leak detector |
| Electrical | insulation tester (megger 500/1000 V), clamp meter, multimeter, phase-rotation tester |
| Pressure test | hydrostatic test pump: water piping at 1.5 × working pressure (≥ 10 bar), hold 2–24 h; refrigerant nitrogen test at design pressure (e.g. ≈ 40 bar for R410A/R32), hold 24 h |

### 22.2 Testing, adjusting & balancing (TAB) and commissioning

| Instrument | Measures | Typical accuracy |
|---|---|---|
| Capture hood (balometer) | diffuser/grille airflow | ±3 % + 7 L/s |
| Thermo-anemometer (vane / hot-wire) | air velocity, temperature | ±2–3 % |
| Pitot tube + digital manometer | duct velocity pressure (airflow traverse), static pressure | ±1 % FS |
| Balancing-valve differential pressure meter | water flow through DRV/PICV | ±1 % |
| Ultrasonic clamp-on flow meter | pipe water flow (non-intrusive) | ±1–2 % |
| Psychrometer / RH-T meter | DB, WB, RH | ±0.2 K, ±2 % RH |
| Contact / IR thermometer, thermal camera | surface and pipe temperatures, insulation defects | ±0.5–1 K |
| CO₂ meter | ventilation verification | ±50 ppm |
| Sound level meter (Class 1/2, octave bands) | NC/NR verification | ±1 dB |
| Optical tachometer | fan/pump speed | ±0.05 % |
| Vibration meter | rotating equipment (ISO 10816) | — |
| Power analyser | kW, PF, harmonics (VSDs) | ±1 % |
| Data loggers (T/RH/kW) | 1–2 week trends, functional performance tests | — |

**TAB sequence:** pre-checks (filters clean, dampers open, fan rotation correct) → main airflow by pitot
traverse → proportional balancing of branches and terminals (target ±10 %) → water-side balancing
(index circuit first, proportional method) → record ESP, motor current and speeds → controls
functional tests (every sequence of operation) → seasonal tests → as-built report.

---

# PART C — EQUIPMENT ENGINEERING & MANUFACTURING

Part C covers how the hardware itself is designed, built and verified. Use it to design or fabricate
components, to check a manufacturer's selection, or to write a technical specification.

## C1. Finned-tube cooling coil: design, rating and fabrication

### C1.1 Standard geometries (plate fin, staggered round tube)

| Tube OD | Tube pitch P_t × row pitch P_l | Tube wall (Cu) | Typical use |
|---|---|---|---|
| 7 mm | 21 × 12.7 mm | 0.25–0.30 mm (inner-grooved) | DX indoor units, VRF |
| 3/8" (9.52 mm) | 25.4 × 22.0 mm | 0.30–0.41 mm | FCU, AHU CHW/DX coils |
| 1/2" (12.7 mm) | 31.75 × 27.5 mm | 0.41–0.51 mm | AHU CHW / HW coils |
| 5/8" (15.88 mm) | 38.1 × 33.0 mm | 0.51–0.89 mm | large AHU, steam coils |

| Fin | Data |
|---|---|
| Material | aluminium 1100/8011 (k ≈ 204–222 W/mK); copper (386 W/mK) for corrosive or high-duty use |
| Thickness | 0.10–0.15 mm (Al), 0.10–0.13 mm (Cu) |
| Fin density | 8–14 FPI (1.8–3.2 mm pitch); keep ≤ 12 FPI for wet coils (carry-over, fouling, cleaning) |
| Pattern | plain; wavy (+20–40 % h); louvered/slit (+50–80 % h, more ΔP and fouling) |
| Coatings | hydrophilic (wet coils, better drainage); epoxy / polyurethane / "blue fin" for coastal and industrial sites; post-coating (e.g. phenolic) for severe corrosion |

### C1.2 Fabrication process

1. **Fin stamping:** a high-speed press with progressive dies punches collared holes (collar height = fin pitch).
2. **Tube preparation:** hairpin bending of the Cu tubes (inner-grooved for refrigerant), then cutting to length.
3. **Stacking:** fins are stacked on a lacing table and the hairpins inserted.
4. **Mechanical expansion:** a bullet (mandrel) is forced through each tube, expanding it about 0.2–0.4 mm in
   diameter to give a tight metal-to-metal bond with the fin collars. Thermal contact resistance should be
   below 0.00002–0.00005 m²K/W; poor expansion is a major hidden capacity loss.
5. **Return bends and headers:** U-bends brazed with BCuP (Cu–P, 0–15 % Ag) filler, under nitrogen purge
   to avoid internal oxide. Headers are Cu or steel, with vent and drain connections.
6. **Casing:** galvanised steel or SS304 end plates and channels, with tube-sheet collars to protect the tubes from abrasion.
7. **Leak test:** pressure test with dry air or N₂ at 2.0–3.0 MPa under water (bubble test), or a helium
   mass-spectrometer test (≤ 1–3 g/year equivalent) for refrigerant coils.
8. **Cleaning and coating:** degrease, apply the coating, then dry.
9. **Rating:** AHRI 410 / EN 1216 test-room certification, or validated rating software.

### C1.3 Rating theory (as implemented in `rate_coil`)

**Air side, Gray & Webb (1986) plain plate fin**, valid for 500 ≤ Re_Dc ≤ 24 700:
```
G_c = ṁ_a / A_min ,   Re_Dc = G_c D_c / μ        (D_c = tube OD + 2 t_fin, the collar diameter)
j₄  = 0.14 Re_Dc^−0.328 (P_t/P_l)^−0.502 (s/D_c)^0.0312       (N ≥ 4 rows; s = fin spacing)
j_N = j₄ · 0.991 [2.24 Re^−0.092 (N/4)^−0.031]^(0.607(4−N))    (N < 4)
h_o = j G_c c_p / Pr^(2/3)
```
**Fin efficiency, Schmidt equivalent circular fin** (staggered tubes):
```
X_M = P_t/2 ,  X_L = ½ √((P_t/2)² + P_l²)
R_eq/r = 1.27 (X_M/r) √(X_L/X_M − 0.3) ,   φ = (R_eq/r − 1)(1 + 0.35 ln(R_eq/r))
m = √(2h/(k t)) ,   η_f = tanh(m r φ)/(m r φ)
Wet fin: use h·c_s/c_p in m (Threlkeld)
```
**Water side:** Gnielinski (§A4.2), with the water velocity set by the number of circuits.

**Wet coil, Braun–Klein–Mitchell enthalpy-effectiveness model.** The driving potential is air enthalpy minus
the enthalpy of saturated air at the water temperature:
```
c_s  = [h_s(t_w,out) − h_s(t_w,in)] / (t_w,out − t_w,in)        (iterate on t_w,out)
R*   = c_s/(h_i A_i) + c_s R″_f/A_i + c_p,a/(h_o η_o,wet A_o)
C_a  = ṁ_a ,  C_w = ṁ_w c_p,w / c_s ,  NTU* = 1/(R* C_min) ,  ε from counter-flow formula
Q    = ε C_min [h_a,in − h_s(t_w,in)]
Leaving air: NTU_o = h_o η_o A_o/(ṁ_a c_p,a),  h_s,eff = h_a,in − (h_a,in − h_a,out)/(1 − e^(−NTU_o))
             t_a,out = t_s,eff + (t_a,in − t_s,eff) e^(−NTU_o),  W_out from h_a,out and t_a,out
```
If the effective surface temperature is above the entering dew point, the coil is dry. The engine then
switches to sensible ε-NTU with UA.

**Water pressure drop:** `Δp = [f L_circ/D_i + K_bend (n_bends)] ½ρV²`, with K_bend ≈ 1.5 for close-radius return bends.
**Air pressure drop:** the engine uses an empirical estimate (±30 %),
`Δp ≈ N_rows · 3.6 · (FPI/12)^0.5 · V_face^1.75 × 1.4 (wet)`.
For accurate air ΔP use the Wang–Chi–Chang (2000) friction correlation or the manufacturer's software.

### C1.4 Design procedure (engine `design_coil`)

1. **Duty:** airflow, entering DB/WB, required total kW and leaving DB (from §9), and the CHW temperatures.
2. **Face area:** `A_face = V̇ / 2.5 m/s`. Fit tubes per row as an integer: `height = n_t × P_t`.
3. **Water flow:** `Q/(ρ c_p ΔT_w)`. Number of circuits for a tube velocity of 0.6–1.5 m/s
   (below about 0.6 m/s flow goes laminar and h_i collapses; above about 1.8 m/s erosion and ΔP rise).
4. **Rate rows × FPI combinations.** Pick the cheapest one (fewest rows × FPI) that meets capacity *and*
   leaving DB, with water ΔP ≤ 40–60 kPa and air ΔP ≤ 150–200 Pa.
5. **Check:** condensate carry-over (face velocity ≤ 2.5 m/s, or add an eliminator), drain-pan size,
   and the laminar-flow risk at part load.

**Worked example** (entering 25.24 °C / 10.47 g/kg, 812 L/s, CHW 6 → 12 °C, 0.69 L/s, 900 × 381 mm face,
3/8" tubes, 15 tubes/row, 11 circuits, Al fins 0.115 mm), from `rate_coil`:

| Rows | FPI | Q kW | SHR | LAT °C | W_out g/kg | t_w,out °C | Δp_w kPa | Δp_air Pa |
|---|---|---|---|---|---|---|---|---|
| 4 | 12 | 13.38 | 0.71 | 15.5 | 8.83 | 10.6 | 13.8 | 91 |
| 6 | 10 | 15.73 | 0.69 | 13.9 | 8.46 | 11.4 | 21.1 | 125 |
| 6 | 12 | 16.63 | 0.70 | 13.2 | 8.36 | 11.7 | 21.1 | 137 |
| 6 | 14 | 17.37 | 0.70 | 12.7 | 8.28 | 12.0 | 21.1 | 148 |
| **8** | **10** | **18.05** | 0.69 | **12.4** | 8.11 | 12.2 | 28.4 | 166 |
| 8 | 12 | 18.93 | 0.69 | 11.8 | 8.00 | 12.5 | 28.4 | 182 |

**Required: 17.46 kW, leaving air ≤ 12.5 °C, so select 8 rows × 10 FPI.**
Air side: Re = 2521, h_o = 56.5 W/m²K, wet fin efficiency 0.69. Water side: 1.03 m/s, Re = 6820,
h_i = 4084 W/m²K. A_o/A_i = 14.8. The air side controls the design; this is why fins exist.
Design levers: CHW at 5 °C, a lower face velocity, or louvered fins would make 6 rows enough.

## C2. Fans: engineering & construction

| Item | Engineering data |
|---|---|
| Theory | Euler: Δp_t = ρ u₂ c_u2 (no inlet swirl). Backward-curved blades (β₂ < 90°) give a stable, non-overloading power curve |
| Efficiency standards | ISO 12759 FMEG; EU Reg. 327/2011 minimum efficiency grades; motor IEC 60034-30-1 (IE3/IE4; EC motors ≈ IE4–IE5) |
| Testing | ISO 5801 / AMCA 210 (performance); ISO 5136 / AMCA 300 (sound power) |
| Balancing | ISO 21940-11 grade G6.3 for general fans, G2.5 for high-speed or critical fans |
| Vibration | ISO 14694 category BV-3 (HVAC): about 4.5 mm/s rms at installation, alarm at 7.1, shutdown at 9.0 (rigid mount) |
| Bearings | L10 life ≥ 40 000 h (ISO 281), ≥ 100 000 h for critical duty |
| Drives | direct or EC (no loss); V-belt 3–5 % loss (flat or synchronous belt 1–2 %) |
| Selection | operate at 80–110 % of peak-efficiency flow, right of the stall line; check parallel-fan instability |

## C3. Pumps: engineering & construction

| Item | Engineering data |
|---|---|
| Euler head | `H = (u₂ c_u2 − u₁ c_u1)/g` |
| NPSH available | `NPSH_a = (p_s − p_v)/(ρg) + z_s − h_f,suction` |
| NPSH margin | NPSH_a ≥ NPSH_r(3 %) + 0.5–1.5 m, or ≥ 1.1–1.3 × NPSH_r |
| Example (engine) | open tower basin 1.5 m above the pump, 32 °C water, 0.8 m suction loss: NPSH_a = (101.3 − 4.76) kPa/(ρg) + 1.5 − 0.8 = **10.6 m** |
| Operating region | preferred 70–120 % of best-efficiency flow (Q_BEP); minimum continuous flow about 25–30 % of Q_BEP |
| Parallel pumps | identical models only; each pump must stay on its curve at single-pump run-out |
| Construction | end-suction (ISO 2858 / 5199), split-case (> 100 L/s), in-line (small); cast-iron casing, bronze or SS impeller, mechanical seal (carbon/SiC), PN16 |
| Testing | ISO 9906 (grade 1B acceptance for critical duty, 2B standard) |
| Curve model | `H = H₀ − k Q²`, fitted from 3 catalogue points; the VSD curve follows the affinity laws |

## C4. Ductwork construction & leakage

| Item | Specification |
|---|---|
| Material | galvanised steel EN 10346 DX51D + Z275 / ASTM A653 G90 (275 g/m² zinc, both sides combined); SS316 for kitchen exhaust or wet areas; aluminium for pools |
| Pressure classes (SMACNA) | ±125, 250, 500, 750, 1000, 1500, 2500 Pa (0.5–10 in.wg) |
| Longitudinal seams | Pittsburgh lock (all pressures up to 1000 Pa); button-punch snap lock (≤ 500 Pa); grooved seam (round) |
| Transverse joints | S-slip and drive (small ducts, ≤ 500 Pa); TDC/TDF roll-formed flange with corner pieces and cleats (≤ 1000–1500 Pa); angle-iron flange (high pressure) |
| Reinforcement | panel width, thickness and pressure class set the intervals (SMACNA tables); cross-break or bead unreinforced panels wider than about 480 mm |
| Sealing | SMACNA seal class A (all seams + joints + penetrations), B (seams + joints), C (joints) |
| Developed area | rectangular: `2(a+b)L` + seam allowance (Pittsburgh ≈ 25–30 mm per seam, TDC ≈ 35–40 mm per flange); round: `πDL`; add 10–15 % for fittings, scrap and laps |
| Hangers | threaded rod + angle trapeze or strap, at ≤ 2.4–3.0 m spacing; supports within 600 mm of each fitting |
| Fire dampers | at fire-compartment walls and floors: EN 15650 / UL 555, rated 1–2 h, access door adjacent |
| Flexible duct | UL 181 / EN 13180; ≤ 1.5 m, no sharp bends |
| Tolerances | ±3 mm on cross-section; flange flatness ±2 mm; diagonal difference ≤ 3 mm |

**Air-tightness, EN 12237 (round) / EN 1507 (rectangular):** `f_max = K · p^0.65` [L/(s·m²)]

| Class | A | B | C | D |
|---|---|---|---|---|
| K | 0.027 | 0.009 | 0.003 | 0.001 |
| Use | not recommended | standard comfort | high-performance / VAV / long runs | cleanrooms, labs |

Example (engine `duct_leakage_en12237`): 42.9 m² of duct at 400 Pa, class B, may leak at most **18.95 L/s**
(2.3 % of 812 L/s). Test sections of ≥ 10 % of the installed area with a fan, a flow-measuring orifice and a
manometer at the design pressure. SMACNA equivalent: leakage class `CL = Q/ΔP^0.65` (cfm/100 ft² at in.wg).
Sealed rectangular ductwork achieves about CL 6, round about CL 3.

## C5. Air-handling-unit construction (EN 1886 / EN 13053)

| Property (model box test) | Class limits |
|---|---|
| Casing strength (deflection at ±1000 Pa) | D1 ≤ 4 mm/m · D2 ≤ 10 mm/m · D3 > 10 mm/m |
| Casing air leakage at −400 Pa | L1 ≤ 0.15 · L2 ≤ 0.44 · L3 ≤ 1.32 L/(s·m²) |
| Thermal transmittance U | T1 ≤ 0.5 · T2 ≤ 1.0 · T3 ≤ 1.4 · T4 ≤ 2.0 W/(m²K) · T5 no requirement |
| Thermal bridging factor k_b | TB1 0.75–1.0 · TB2 0.6–0.75 · TB3 0.45–0.6 · TB4 0.3–0.45 · TB5 none |
| Filter bypass leakage | ≤ 0.5 % (F9) up to 6 % (M5) of flow |

Specify **D1 / L1 / T2 / TB2** for humid climates, which avoids external condensation on casings.
Construction details:
- **Panels:** double-skin 50 mm panels, PU (k ≈ 0.022) or mineral wool (k ≈ 0.036, non-combustible).
- **Frame:** thermally broken aluminium.
- **Drain pan:** SS304, double-sloped about 1:100, connection at the lowest point, extending past the coil plus eliminator.
- **Coil handling:** coils on rails for withdrawal.
- **Fans:** plug or EC fans on anti-vibration mounts, with a hinged access door and window.
- **Other:** internal light, and pressure test points across filters and coils.

| Heat recovery | Sensible ε | Latent | Notes |
|---|---|---|---|
| Enthalpy wheel | 0.70–0.80 | yes (0.60–0.75) | 1–5 % carry-over; purge sector |
| Plate, cross-flow | 0.50–0.70 | no (membrane type: yes) | no cross-contamination |
| Plate, counter-flow | 0.75–0.90 | no | frost control in cold climates |
| Run-around coil | 0.45–0.55 | no | separated airstreams |
| Heat pipe | 0.40–0.60 | no | wrap-around for dehumidification + free reheat |

**Filtration classes** (ISO 16890 vs EN 779 vs ASHRAE 52.2; the mapping is approximate):

| ISO 16890 | EN 779 (old) | MERV | Use |
|---|---|---|---|
| ISO Coarse 60–80 % | G3–G4 | 6–8 | pre-filter |
| ePM10 50–65 % | M5 | 9–10 | general |
| ePM2.5 50–65 % | M6 | 11–12 | offices (min. supply) |
| ePM1 50–70 % | F7 | 13–14 | offices, hotels (final) |
| ePM1 80–95 % | F8–F9 | 15–16 | hospitals, urban pollution |
| HEPA H13 / H14 (EN 1822) | — | — | 99.95 / 99.995 % at MPPS: cleanrooms, operating theatres |

## C6. Thermal insulation engineering

**Condensation criterion:** the outer jacket surface must stay above the ambient dew point plus a margin.
```
q = (t_amb − t_fluid) / (R_ins + R_s) ,  R_ins = ln(r₂/r₁)/(2πk) ,  R_s = 1/(h_s 2π r₂)
t_surface = t_amb − q·R_s ≥ t_dp + 0.5 K       → solve for r₂ (engine insulation_condensation)
```
h_s ≈ 9 W/m²K for high-emissivity surfaces (elastomeric, painted) and ≈ 5.7 W/m²K for bright aluminium or
SS jackets (low ε, so less radiation). **Low-emissivity jackets need thicker insulation.**

Engine result for CHW at 6 °C, ambient 32 °C / 80 % RH (dew point 28.1 °C), k = 0.036 W/mK:

| Pipe | OD mm | min. thickness, h_s = 9 | min. thickness, h_s = 5.7 (Al jacket) | Specify |
|---|---|---|---|---|
| DN25 | 33.4 | 18.8 mm | 27.1 mm | 25 mm (32 with jacket) |
| DN50 | 60.3 | 20.8 mm | 30.4 mm | 25 mm (32 with jacket) |
| DN100 | 114.3 | 22.7 mm | 33.7 mm | 32 mm (40 with jacket) |
| DN150 | 168.3 | 23.6 mm | 35.5 mm | 32 mm (40 with jacket) |
| DN200 | 219.1 | 24.2 mm | 36.6 mm | 40 mm (40–50 with jacket) |

The vapour barrier is as important as the thickness. A single unsealed joint lets moisture migrate in,
the insulation becomes wet, k rises 2–5×, and pipes corrode. Elastomeric foam needs μ ≥ 7000; glass wool
always needs a sealed foil (FSK) vapour retarder.

| Material | k at 10 °C [W/mK] | Vapour resistance μ | Service temp. °C | Fire |
|---|---|---|---|---|
| Elastomeric (NBR/EPDM) | 0.033–0.036 | 7 000–10 000 | −50 … +105 | B-s3,d0 / BL-s3,d0 |
| PIR / PUR | 0.022–0.026 | 50–100 (needs barrier) | −180 … +120 | B–E |
| Phenolic | 0.018–0.021 | 35–50 | −180 … +120 | B-s2,d0 |
| Glass wool | 0.032–0.040 | ≈ 1 (needs FSK) | … +250 | A1/A2 |
| Rock wool | 0.034–0.040 | ≈ 1 (needs barrier) | … +650 | A1 |
| Cellular glass | 0.038–0.045 | ∞ | −260 … +430 | A1 |

## C7. Acoustics

```
Level addition:  L = 10 log₁₀ Σ 10^(Lᵢ/10)        (2 equal sources: +3 dB; 10 equal: +10 dB)
Room equation:   L_p = L_w + 10 log₁₀[ Q/(4π r²) + 4/R ] ,  R = S ᾱ/(1 − ᾱ)   (room constant, m²)
Directivity Q:   1 free field, 2 on a floor/wall, 4 at an edge, 8 in a corner
```
Engine example: sources of 75 + 70 + 70 dB give 77.1 dB Lw. At 3 m in a room with R = 200 m²,
**L_p = 62.9 dB**, so the equipment needs a plant room or silencers.

| Path element | Attenuation |
|---|---|
| Unlined sheet-metal duct | 0.1–0.6 dB/m (low frequencies best) |
| Lined duct, 25–50 mm | 1–10 dB/m at mid frequencies |
| Lined elbow | 2–10 dB |
| Duct end reflection (small outlets) | 5–15 dB at 63–250 Hz |
| Rectangular / round silencer, 1–2 m | insertion loss 10–35 dB (mid frequencies), self-noise at > 10 m/s face velocity |

Regenerated noise rises with about 50–60 log₁₀(V), so keep velocities to §12.3. Check break-out noise
from rectangular ducts above quiet rooms.

## C8. Vibration isolation

Single-degree-of-freedom system: `f_n = (1/2π)√(k/m) = 15.76/√δ_st[mm]` Hz; transmissibility (undamped)
`T = 1/|r² − 1|`, with `r = f_forcing/f_n`. Isolation starts at r > √2; aim for r ≥ 3–5 (T ≤ 0.12–0.04).
Engine example: a 1450 rpm fan (24.2 Hz) on 25 mm springs has f_n = 3.15 Hz, r = 7.7, T = 0.017,
so **98.3 % isolation**.

| Equipment / speed | Minimum static deflection (ground floor → upper floor, long span) |
|---|---|
| > 1500 rpm | 6–10 mm → 25 mm |
| 900–1500 rpm | 10–25 mm → 25–50 mm |
| 500–900 rpm (towers, large fans) | 25 mm → 50–75 mm |
| Chillers (screw/centrifugal) | 6 mm neoprene → 25–50 mm springs |
| Pumps | inertia base (1.5–2.5 × pump mass) + 25 mm springs, flexible connectors |

Isolator types: neoprene pads (2–5 mm), neoprene mounts (5–10 mm), open springs (25–75 mm), restrained
springs (outdoor, seismic), spring hangers for suspended units and pipes near equipment.

## C9. DX / VRF refrigerant installation practice

- **Pipe material:** refrigeration-grade copper (EN 12735-1 / ASTM B280), cleaned and capped, with wall
  thickness for the design pressure. R410A and R32 run at about 4.2 MPa, so check wall ratings, especially ≥ 3/4" annealed tube.
- **Brazing:** braze with nitrogen flowing at 2–5 L/min (prevents oxide scale that blocks EEVs and filters).
  Flare joints only where the manufacturer permits.
- **Strength / tightness test:** N₂ at the manufacturer's design pressure (e.g. 4.15 MPa), hold 24 h;
  correct for temperature (≈ 0.01 MPa per K).
- **Evacuation:** triple evacuation or a deep vacuum to ≤ 500 micron (67 Pa). Isolate and hold 30–60 min;
  a rise to ≤ 1000 micron is acceptable, a continuous rise means a leak, and a plateau means moisture.
- **Charge:** factory charge + additional charge per metre of liquid line (manufacturer table), weighed on a scale.
- **Refrigerant concentration limit (EN 378 / ISO 5149):** charge/room volume ≤ practical limit:
  R410A 0.44 kg/m³, R134a 0.25 kg/m³, R32 0.061 kg/m³ (A2L, flammability-based). Otherwise add leak
  detection + ventilation or shut-off valves, or redesign the zoning.
- **Oil return:** follow §A6.3 velocities; traps at the base of suction risers and about every 7.5 m of rise
  for non-VRF systems. VRF systems follow the manufacturer's height and length limits.

---

# PART D — STANDARDS, QUALITY ASSURANCE, OPERATION

## D1. Standards map

| Area | Standards |
|---|---|
| Design, loads & comfort | ASHRAE Handbook Fundamentals (Ch.1, 14, 18), ASHRAE 55, ISO 7730, CIBSE Guide A, EN 16798-1 |
| Ventilation & IAQ | ASHRAE 62.1 / 62.2, EN 16798-3, ASHRAE 170 (healthcare), local codes |
| Energy | ASHRAE 90.1, EN 16798, national energy codes (e.g. QCVN 09:2017/BXD) |
| Ductwork | SMACNA HVAC Duct Construction Standards, EN 1505/1506/1507, EN 12237, DW/144 (UK) |
| Piping | ASME B31.9 (building services), ASME B31.5 (refrigeration), EN 13480 |
| Refrigerant safety | ASHRAE 15 & 34, ISO 5149, ISO 817, EN 378 |
| Equipment rating | AHRI 550/590 (chillers), 410 (coils), 440 (FCU), 210/240 (unitary), 1230 (VRF), EN 14511, Eurovent certification |
| AHU | EN 1886, EN 13053, ISO 16890 (filters), EN 1822 (HEPA), ASHRAE 52.2 |
| Fans & pumps | ISO 5801, AMCA 210/300, ISO 12759, ISO 9906, ISO 5199, IEC 60034-30-1 |
| Fire & smoke | NFPA 90A/92, EN 12101, UL 555 / EN 15650 (fire dampers), local fire code |
| TAB & commissioning | ASHRAE 111, NEBB, AABC, CIBSE Code A/W/C, ASHRAE Guideline 0 / Standard 202 |
| Legionella | ASHRAE 188, CTI guidelines, local health rules |

## D2. Inspection & test plan (ITP), with acceptance criteria

| Stage | Test | Acceptance |
|---|---|---|
| Factory | chiller FAT (AHRI 550/590) | capacity and kW/TR within the standard tolerance (≈ ±5 % at full load, tolerance formula) |
| Factory | AHU casing / fan run test | EN 1886 classes as specified; vibration ≤ ISO 14694 BV-3 |
| Factory | coil leak test | no bubbles at 2.0–3.0 MPa / helium criterion |
| Site | duct leakage test (≥ 10 % of area) | ≤ EN 12237 class (e.g. B) or SMACNA CL |
| Site | hydrostatic pipe test | 1.5 × working pressure (≥ 10 bar), 2 h (or per spec), no pressure drop after temperature correction, no visible leaks |
| Site | refrigerant N₂ test and vacuum | 24 h at design pressure with no drop; vacuum ≤ 500 micron, hold test passed |
| Site | flushing & cleaning (BSRIA BG29) | filter mesh clean, water TSS / iron within spec before connecting coils |
| Site | insulation inspection | thickness as specified, vapour barrier continuous, no gaps at hangers (use load-bearing inserts) |
| Site | electrical | insulation resistance ≥ 1 MΩ (megger), rotation, overload settings, earthing |
| TAB | airflow | terminals ±10 %, main ducts ±5–10 % of design |
| TAB | water flow | coils ±10 %, plant ±5 % |
| Performance | room conditions | 24 ±1 °C, RH ≤ 60 %, NC target, CO₂ ≤ outdoor + 700 ppm |
| Controls | point-to-point and functional tests | every BMS point verified; every sequence demonstrated (setpoint, alarms, failure modes) |
| Handover | documents | as-built drawings, O&M manuals, TAB report, Cx report, training records, warranties |

## D3. Operation & maintenance (preventive schedule)

| Item | Frequency | Action / trigger |
|---|---|---|
| Filters | by ΔP (change at final ΔP) or 1–3 months | record ΔP weekly on the BMS |
| Coils & drain pans | 6–12 months | clean fins (low-pressure water/foam), check drainage and traps; biocide tablets in humid climates |
| Fans & belts | 3–6 months | belt tension and alignment, bearing lubrication, vibration trend |
| Pumps | 6 months | seal leakage, vibration, motor current vs. curve |
| Chillers | daily log; annual service | approach temperatures: evaporator > 2 K or condenser > 2–3 K above baseline means tube cleaning; oil analysis; refrigerant leak check |
| Cooling towers | weekly | water treatment (conductivity/TDS, pH, biocide), drift eliminators; **Legionella** risk plan, culture tests monthly to quarterly, clean and disinfect at least 2×/year |
| Water treatment, closed loop | quarterly | inhibitor level, pH 8.5–9.5, iron, bacteria |
| Controls & sensors | annually | calibrate T/RH/CO₂/ΔP sensors; review trends (simultaneous heating and cooling, hunting loops, overrides) |
| Energy | monthly | kW/TR of the plant, SFP of fans, EUI benchmark |
