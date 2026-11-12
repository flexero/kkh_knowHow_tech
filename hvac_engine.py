#!/usr/bin/env python3
"""
hvac_engine.py - HVAC design calculation engine (SI units, Python stdlib only)

Modules
  1. Psychrometrics          (ASHRAE Fundamentals 2017/2021, Ch.1)
  2. Climate & solar         (ASHRAE clear-sky model, daily temperature profile)
  3. Cooling load            (24 h hourly: sol-air + decrement/time-lag conduction,
                              solar through glass, internal gains, radiant time delay)
  4. Heating load            (steady state)
  5. Ventilation             (ASHRAE 62.1 Ventilation Rate Procedure)
  6. Air-side process        (supply airflow, mixing, coil ADP / bypass factor)
  7. Duct sizing             (Darcy-Weisbach + Colebrook, equal friction)
  8. Pipe sizing             (chilled / condenser water)
  9. Plant & power           (chiller, cooling tower, fan & pump power, fan laws)

Units: temperature degC, pressure kPa (psychrometrics) or Pa (duct/pipe),
humidity ratio kg/kg dry air, enthalpy kJ/kg dry air, power W, airflow m3/s.

Run `python hvac_engine.py` to execute the worked example (office zone).
"""

# ---------------------------------------------------------------------------
# Compatibility block: runs on Python 3.6 ... 3.15 with NO third-party packages.
# (f-strings require >= 3.6; dataclasses are emulated below on 3.6.)
# ---------------------------------------------------------------------------
import sys
if sys.version_info < (3, 6):
    sys.exit("This engine needs Python 3.6 or newer (found %d.%d)." % sys.version_info[:2])
from typing import Dict, List, Optional, Tuple  # noqa: F401  (typing exists since 3.5)
try:
    from dataclasses import dataclass, field     # Python >= 3.7
except ImportError:                              # Python 3.6: minimal built-in fallback
    _MISSING = object()

    class _Field(object):
        def __init__(self, default=_MISSING, default_factory=_MISSING):
            self.default, self.default_factory = default, default_factory

    def field(default=_MISSING, default_factory=_MISSING):
        return _Field(default, default_factory)

    def dataclass(cls):
        names = list(cls.__dict__.get("__annotations__", {}))
        defaults = {n: cls.__dict__[n] for n in names if n in cls.__dict__}

        def __init__(self, *args, **kw):
            if len(args) > len(names):
                raise TypeError("%s: too many arguments" % cls.__name__)
            for i, n in enumerate(names):
                if i < len(args):
                    v = args[i]
                elif n in kw:
                    v = kw.pop(n)
                elif n in defaults:
                    d = defaults[n]
                    if isinstance(d, _Field):
                        v = d.default_factory() if d.default_factory is not _MISSING else d.default
                    else:
                        v = d
                else:
                    raise TypeError("%s: missing argument '%s'" % (cls.__name__, n))
                setattr(self, n, v)
            if kw:
                raise TypeError("%s: unexpected argument(s) %s" % (cls.__name__, list(kw)))

        def __repr__(self):
            return "%s(%s)" % (cls.__name__, ", ".join("%s=%r" % (n, getattr(self, n)) for n in names))

        def __eq__(self, other):
            return type(other) is type(self) and all(getattr(self, n) == getattr(other, n) for n in names)

        for n, d in defaults.items():           # remove Field placeholders from the class
            if isinstance(d, _Field):
                if d.default is not _MISSING:
                    setattr(cls, n, d.default)
                else:
                    delattr(cls, n)
        cls.__init__, cls.__repr__, cls.__eq__ = __init__, __repr__, __eq__
        return cls
# ---------------------------------------------------------------------------
import math

# ============================================================================
# 1. PSYCHROMETRICS
# ============================================================================
R_DA = 0.287042          # kJ/(kg.K) gas constant dry air
MW_RATIO = 0.621945      # Mw/Mda
CP_DA = 1.006            # kJ/(kg.K) dry air
CP_V = 1.86              # kJ/(kg.K) water vapour
HFG0 = 2501.0            # kJ/kg latent heat at 0 degC


def pressure_at_altitude(z_m: float) -> float:
    """Standard atmospheric pressure [kPa] at altitude z [m] (ASHRAE eq.3)."""
    return 101.325 * (1.0 - 2.25577e-5 * z_m) ** 5.2559


def p_ws(t: float) -> float:
    """Saturation pressure of water vapour [kPa] at t [degC] (Hyland-Wexler)."""
    T = t + 273.15
    if t >= 0.0:  # over liquid water, 0..200 degC
        ln = (-5.8002206e3 / T + 1.3914993 - 4.8640239e-2 * T
              + 4.1764768e-5 * T ** 2 - 1.4452093e-8 * T ** 3 + 6.5459673 * math.log(T))
    else:         # over ice, -100..0 degC
        ln = (-5.6745359e3 / T + 6.3925247 - 9.6778430e-3 * T + 6.2215701e-7 * T ** 2
              + 2.0747825e-9 * T ** 3 - 9.4840240e-13 * T ** 4 + 4.1635019 * math.log(T))
    return math.exp(ln) / 1000.0


def w_from_pw(pw: float, p: float = 101.325) -> float:
    """Humidity ratio [kg/kg] from vapour partial pressure [kPa]."""
    return MW_RATIO * pw / (p - pw)


def pw_from_w(w: float, p: float = 101.325) -> float:
    return p * w / (MW_RATIO + w)


def w_sat(t: float, p: float = 101.325) -> float:
    return w_from_pw(p_ws(t), p)


def w_from_rh(t: float, rh: float, p: float = 101.325) -> float:
    """rh as fraction 0..1."""
    return w_from_pw(rh * p_ws(t), p)


def w_from_wb(t: float, twb: float, p: float = 101.325) -> float:
    """Humidity ratio from dry bulb and thermodynamic wet bulb (ASHRAE eq.33/35)."""
    ws = w_sat(twb, p)
    if twb >= 0:
        return ((2501 - 2.326 * twb) * ws - 1.006 * (t - twb)) / (2501 + 1.86 * t - 4.186 * twb)
    return ((2830 - 0.24 * twb) * ws - 1.006 * (t - twb)) / (2830 + 1.86 * t - 2.1 * twb)


def enthalpy(t: float, w: float) -> float:
    """Moist air enthalpy [kJ/kg dry air]."""
    return CP_DA * t + w * (HFG0 + CP_V * t)


def t_from_hw(h: float, w: float) -> float:
    return (h - HFG0 * w) / (CP_DA + CP_V * w)


def specific_volume(t: float, w: float, p: float = 101.325) -> float:
    """[m3/kg dry air]"""
    return R_DA * (t + 273.15) * (1 + 1.607858 * w) / p


def rh_from_w(t: float, w: float, p: float = 101.325) -> float:
    return pw_from_w(w, p) / p_ws(t)


def _bisect(f, lo: float, hi: float, tol: float = 1e-6, it: int = 200) -> float:
    flo = f(lo)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < tol or (hi - lo) < tol:
            return mid
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def dew_point(w: float, p: float = 101.325) -> float:
    pw = pw_from_w(w, p)
    return _bisect(lambda t: p_ws(t) - pw, -60.0, 100.0)


def wet_bulb(t: float, w: float, p: float = 101.325) -> float:
    return _bisect(lambda twb: w_from_wb(t, twb, p) - w, -60.0, t)


@dataclass
class AirState:
    """Complete moist-air state. Build with AirState.from_*()."""
    t: float      # dry bulb degC
    w: float      # humidity ratio kg/kg
    p: float = 101.325

    @classmethod
    def from_rh(cls, t, rh, p=101.325):
        return cls(t, w_from_rh(t, rh, p), p)

    @classmethod
    def from_wb(cls, t, twb, p=101.325):
        return cls(t, w_from_wb(t, twb, p), p)

    @property
    def h(self): return enthalpy(self.t, self.w)
    @property
    def rh(self): return rh_from_w(self.t, self.w, self.p)
    @property
    def v(self): return specific_volume(self.t, self.w, self.p)
    @property
    def rho(self): return (1 + self.w) / self.v          # moist air density kg/m3
    @property
    def tdp(self): return dew_point(self.w, self.p)
    @property
    def twb(self): return wet_bulb(self.t, self.w, self.p)

    def __str__(self):
        return (f"DB {self.t:6.2f} C | WB {self.twb:6.2f} C | DP {self.tdp:6.2f} C | "
                f"RH {self.rh*100:5.1f} % | W {self.w*1000:6.2f} g/kg | "
                f"h {self.h:6.2f} kJ/kg | v {self.v:.4f} m3/kg")


def mix(a: AirState, ma: float, b: AirState, mb: float) -> AirState:
    """Adiabatic mixing of two streams by dry-air mass (or volume, approx.)."""
    w = (ma * a.w + mb * b.w) / (ma + mb)
    h = (ma * a.h + mb * b.h) / (ma + mb)
    return AirState(t_from_hw(h, w), w, a.p)


# ============================================================================
# 2. CLIMATE & SOLAR
# ============================================================================
# ASHRAE fraction of daily range, hour 1..24 (solar time). T = Tmax - f*DR
DAILY_RANGE_FRAC = [0.87, 0.92, 0.96, 0.99, 1.00, 0.98, 0.93, 0.84, 0.71, 0.56, 0.39, 0.23,
                    0.11, 0.03, 0.00, 0.03, 0.10, 0.21, 0.34, 0.47, 0.58, 0.68, 0.76, 0.82]


def outdoor_profile(t_db_max: float, daily_range: float, t_wb_mc: float, p=101.325):
    """24 hourly outdoor AirStates. Humidity ratio held constant at design (MCWB)."""
    w = w_from_wb(t_db_max, t_wb_mc, p)
    return [AirState(t_db_max - f * daily_range, w, p) for f in DAILY_RANGE_FRAC]


@dataclass
class Site:
    latitude: float          # deg, north positive
    longitude: float         # deg, east positive
    tz: float                # standard time zone, h from UTC (e.g. +7)
    altitude: float = 0.0    # m
    tau_b: float = 0.45      # ASHRAE clear-sky beam optical depth (month/site)
    tau_d: float = 2.20      # diffuse optical depth
    ground_refl: float = 0.2


def solar_position(site: Site, day_of_year: int, hour_solar: float):
    """Return (altitude beta deg, azimuth phi deg from south, +west) at solar time."""
    n = day_of_year
    decl = 23.45 * math.sin(math.radians(360 * (284 + n) / 365))
    H = 15.0 * (hour_solar - 12.0)                    # hour angle, deg
    L, d, Hr = map(math.radians, (site.latitude, decl, H))
    sb = math.cos(L) * math.cos(d) * math.cos(Hr) + math.sin(L) * math.sin(d)
    beta = math.asin(max(-1, min(1, sb)))
    cb = math.cos(beta)
    if cb < 1e-9:
        return math.degrees(beta), 0.0
    sphi = math.sin(Hr) * math.cos(d) / cb
    cphi = (math.cos(Hr) * math.cos(d) * math.sin(L) - math.sin(d) * math.cos(L)) / cb
    return math.degrees(beta), math.degrees(math.atan2(sphi, cphi))


def clear_sky(site: Site, day_of_year: int, beta_deg: float):
    """ASHRAE clear-sky model -> (Eb normal beam, Ed horizontal diffuse) W/m2."""
    if beta_deg <= 0:
        return 0.0, 0.0
    E0 = 1367 * (1 + 0.033 * math.cos(math.radians(360 * (day_of_year - 3) / 365)))
    m = 1.0 / (math.sin(math.radians(beta_deg)) + 0.50572 * (6.07995 + beta_deg) ** -1.6364)
    tb, td = site.tau_b, site.tau_d
    ab = 1.454 - 0.406 * tb - 0.268 * td + 0.021 * tb * td
    ad = 0.507 + 0.205 * tb - 0.080 * td - 0.190 * tb * td
    return E0 * math.exp(-tb * m ** ab), E0 * math.exp(-td * m ** ad)


def surface_irradiance(site: Site, day: int, hour: float, tilt: float, surf_az: float):
    """Total irradiance on a surface [W/m2] and its beam part.
    tilt: 90 = wall, 0 = roof.  surf_az: 0=S, 90=W, -90=E, 180=N (deg)."""
    beta, phi = solar_position(site, day, hour)
    Eb, Ed = clear_sky(site, day, beta)
    if beta <= 0:
        return 0.0, 0.0
    S, B, G = map(math.radians, (tilt, beta, phi - surf_az))
    cos_t = math.cos(B) * math.cos(G) * math.sin(S) + math.sin(B) * math.cos(S)
    Et_b = Eb * max(0.0, cos_t)
    Y = max(0.45, 0.55 + 0.437 * cos_t + 0.313 * cos_t ** 2)
    Et_d = Ed * (Y * math.sin(S) + math.cos(S)) if tilt <= 90 else Ed * Y * math.sin(S)
    Et_r = (Eb * math.sin(B) + Ed) * site.ground_refl * (1 - math.cos(S)) / 2
    return Et_b + Et_d + Et_r, Et_b


# ============================================================================
# 3. COOLING LOAD (hourly, 24 h)
# ============================================================================
# Radiant fraction of each heat gain (ASHRAE Fund. Ch.18 Table 14, rounded)
RADIANT_FRAC = {"people": 0.60, "lighting": 0.57, "equipment": 0.20, "wall": 0.46,
                "roof": 0.60, "glass_cond": 0.46, "solar": 0.90}

# Room thermal-storage weight -> geometric decay ratio of the radiant time series
RTS_DECAY = {"light": 0.30, "medium": 0.50, "heavy": 0.70}


def radiant_time_series(weight: str = "medium") -> List[float]:
    """Simplified RTS: r_j = (1-a) a^j normalised over 24 h (sum = 1)."""
    a = RTS_DECAY[weight]
    r = [(1 - a) * a ** j for j in range(24)]
    s = sum(r)
    return [x / s for x in r]


def apply_rts(gain24: List[float], rad_frac: float, rts: List[float]) -> List[float]:
    """Split gain into convective (instant) + radiant (delayed by RTS) -> cooling load."""
    out = []
    for h in range(24):
        rad = sum(rts[j] * gain24[(h - j) % 24] for j in range(24)) * rad_frac
        out.append(gain24[h] * (1 - rad_frac) + rad)
    return out


@dataclass
class Opaque:
    name: str
    area: float           # m2 net
    U: float              # W/(m2.K)
    tilt: float           # 90 wall, 0 roof
    azimuth: float        # 0 S, 90 W, -90 E, 180 N
    absorptance: float = 0.7
    decrement: float = 0.5    # f (0..1)
    lag_h: int = 6            # time lag, h


@dataclass
class Window:
    name: str
    area: float
    U: float
    SHGC: float
    azimuth: float
    tilt: float = 90.0
    shade_coeff: float = 1.0   # internal shading multiplier (IAC), 1 = none


@dataclass
class Zone:
    name: str
    floor_area: float
    height: float
    t_in: float = 24.0
    rh_in: float = 0.50
    walls: list = field(default_factory=list)
    windows: list = field(default_factory=list)
    people: int = 0
    people_sens: float = 75.0       # W/person (moderately active office)
    people_lat: float = 55.0
    lighting_wpm2: float = 9.0
    equipment_wpm2: float = 12.0
    schedule: list = field(default_factory=lambda: [0]*8 + [1]*10 + [0]*6)  # occupied 08-18
    infiltration_ach: float = 0.0
    construction: str = "medium"


def sol_air(t_out: float, I: float, absorp: float, tilt: float, h_o: float = 17.0) -> float:
    """Sol-air temperature [degC]. Long-wave correction 3.9 K for horizontal only."""
    lw = 3.9 if tilt < 45 else 0.0
    return t_out + absorp * I / h_o - lw


def cooling_load(zone: Zone, site: Site, day: int, t_db_max: float, daily_range: float,
                 t_wb_mc: float):
    """Hourly sensible & latent room cooling load [W]. Returns dict of 24-lists."""
    p = pressure_at_altitude(site.altitude)
    out = outdoor_profile(t_db_max, daily_range, t_wb_mc, p)
    room = AirState.from_rh(zone.t_in, zone.rh_in, p)
    rts = radiant_time_series(zone.construction)
    hours = [h + 0.5 for h in range(24)]           # mid-hour solar time
    comp = {}

    # Opaque conduction: q(t) = U A [ (Tsa_mean - Ti) + f (Tsa(t-lag) - Tsa_mean) ]
    wall_tot, roof_tot = [0.0]*24, [0.0]*24
    for s in zone.walls:
        tsa = [sol_air(out[h].t, surface_irradiance(site, day, hours[h], s.tilt, s.azimuth)[0],
                       s.absorptance, s.tilt) for h in range(24)]
        mean = sum(tsa) / 24
        q = [s.U * s.area * ((mean - zone.t_in) + s.decrement * (tsa[(h - s.lag_h) % 24] - mean))
             for h in range(24)]
        tgt = roof_tot if s.tilt < 45 else wall_tot
        for h in range(24):
            tgt[h] += q[h]
    comp["walls"] = apply_rts(wall_tot, RADIANT_FRAC["wall"], rts)
    comp["roof"] = apply_rts(roof_tot, RADIANT_FRAC["roof"], rts)

    # Windows: conduction (instant) + solar SHGC * Et
    gc, gs = [0.0]*24, [0.0]*24
    for wdw in zone.windows:
        for h in range(24):
            gc[h] += wdw.U * wdw.area * (out[h].t - zone.t_in)
            Et, _ = surface_irradiance(site, day, hours[h], wdw.tilt, wdw.azimuth)
            gs[h] += wdw.area * wdw.SHGC * wdw.shade_coeff * Et
    comp["glass_cond"] = apply_rts(gc, RADIANT_FRAC["glass_cond"], rts)
    comp["solar"] = apply_rts(gs, RADIANT_FRAC["solar"], rts)

    # Internal gains with schedule
    sch = zone.schedule
    comp["people"] = apply_rts([zone.people * zone.people_sens * sch[h] for h in range(24)],
                               RADIANT_FRAC["people"], rts)
    comp["lighting"] = apply_rts([zone.lighting_wpm2 * zone.floor_area * sch[h] for h in range(24)],
                                 RADIANT_FRAC["lighting"], rts)
    comp["equipment"] = apply_rts([zone.equipment_wpm2 * zone.floor_area * (0.15 + 0.85*sch[h])
                                   for h in range(24)], RADIANT_FRAC["equipment"], rts)

    # Infiltration (all convective)
    vol = zone.floor_area * zone.height
    q_inf = zone.infiltration_ach * vol / 3600.0       # m3/s
    rho = room.rho
    comp["infiltration_s"] = [q_inf * rho * 1006 * (out[h].t - zone.t_in) for h in range(24)]
    lat_inf = [q_inf * rho * 2501e3 * max(0.0, out[h].w - room.w) for h in range(24)]
    lat_ppl = [zone.people * zone.people_lat * sch[h] for h in range(24)]

    sens = [sum(c[h] for c in comp.values()) for h in range(24)]
    lat = [lat_inf[h] + lat_ppl[h] for h in range(24)]
    peak_h = max(range(24), key=lambda h: sens[h] + lat[h])
    return {"components": comp, "sensible": sens, "latent": lat, "peak_hour": peak_h,
            "outdoor": out, "room": room}


# ============================================================================
# 4. HEATING LOAD (steady state)
# ============================================================================
def heating_load(elements: List[Tuple[float, float]], t_in: float, t_out: float,
                 volume: float = 0.0, ach: float = 0.0, safety: float = 1.15) -> float:
    """elements = [(U, A), ...].  Q = (sum UA + 0.34*ACH*V) * dT * safety  [W]
    0.34 = rho*cp/3600 [Wh/(m3.K)]."""
    ua = sum(U * A for U, A in elements)
    return (ua + 0.34 * ach * volume) * (t_in - t_out) * safety


# ============================================================================
# 5. VENTILATION - ASHRAE 62.1 VRP
# ============================================================================
# Space type: (Rp L/s.person, Ra L/s.m2, default occupant density #/100 m2)
ASHRAE621 = {
    "office":            (2.5, 0.3, 5),
    "conference":        (2.5, 0.3, 50),
    "reception":         (2.5, 0.3, 30),
    "classroom_9plus":   (3.8, 0.6, 35),
    "lecture_hall":      (3.8, 0.3, 65),
    "retail_sales":      (3.8, 0.6, 15),
    "restaurant_dining": (3.8, 0.9, 70),
    "hotel_bedroom":     (2.5, 0.3, 10),
    "health_club":       (10.0, 0.3, 40),
    "lobby":             (2.5, 0.3, 10),
}


def outdoor_air(space: str, area: float, people: Optional[int] = None, Ez: float = 1.0):
    """Return (Vbz, Voz) in L/s.  Vbz = Rp*Pz + Ra*Az ;  Voz = Vbz/Ez."""
    Rp, Ra, dens = ASHRAE621[space]
    Pz = people if people is not None else math.ceil(dens * area / 100)
    Vbz = Rp * Pz + Ra * area
    return Vbz, Vbz / Ez


def system_oa_multizone(zones: List[Tuple[float, float]], Ev_min_note: bool = True):
    """Simplified multi-zone VRP. zones = [(Voz L/s, Vpz primary supply L/s), ...].
    Zpz = Voz/Vpz ; Xs = Vou/Vps ; Ev = 1 + Xs - Zpz_max  (ASHRAE 62.1 App. A simplified)
    Vou approximated as sum(Voz) (no diversity). Returns (Vot, Ev)."""
    Vou = sum(v for v, _ in zones)
    Vps = sum(p for _, p in zones)
    Xs = Vou / Vps
    Zmax = max(v / p for v, p in zones)
    Ev = min(1.0, 1 + Xs - Zmax)
    return Vou / Ev, Ev


# ============================================================================
# 6. AIR-SIDE PROCESS: supply air, mixing, coil
# ============================================================================
def supply_airflow(q_sens: float, t_room: float, t_supply: float, rho: float = 1.2,
                   cp: float = 1026.0) -> float:
    """Supply airflow [m3/s] = Qs / (rho cp dT).  cp 1026 includes vapour."""
    return q_sens / (rho * cp * (t_room - t_supply))


def supply_humidity(room: AirState, q_lat: float, flow: float, rho: float = 1.2) -> float:
    """Required supply humidity ratio to absorb latent load."""
    return room.w - q_lat / (rho * flow * 2501e3)


def coil_adp_bf(entering: AirState, leaving: AirState):
    """Extend coil line to saturation -> apparatus dew point and bypass factor."""
    def f(x):  # x: fraction along line beyond leaving (0..)
        t = leaving.t + x * (leaving.t - entering.t)
        w = leaving.w + x * (leaving.w - entering.w)
        return w - w_sat(t, entering.p)
    x = _bisect(f, 0.0, 5.0)
    adp_t = leaving.t + x * (leaving.t - entering.t)
    bf = (leaving.t - adp_t) / (entering.t - adp_t)
    return adp_t, bf


def coil_load(entering: AirState, leaving: AirState, flow_m3s: float):
    """Coil total, sensible [W] and condensate [kg/s] (flow at entering density)."""
    m = flow_m3s / entering.v                       # kg dry air /s
    total = m * (entering.h - leaving.h) * 1000
    sens = m * (CP_DA + CP_V * entering.w) * (entering.t - leaving.t) * 1000
    cond = m * (entering.w - leaving.w)
    return total, sens, cond


# ============================================================================
# 7. DUCT SIZING
# ============================================================================
ROUGHNESS = {"galv_steel": 0.09e-3, "aluminium": 0.05e-3, "flex_duct": 3.0e-3,
             "fibreglass": 0.9e-3, "steel_pipe": 0.045e-3, "pvc": 0.0015e-3,
             "copper": 0.0015e-3}
ROUND_DUCT_MM = [100, 125, 160, 200, 250, 315, 355, 400, 450, 500, 560, 630, 710,
                 800, 900, 1000, 1120, 1250, 1400, 1600]
AIR_NU = 1.53e-5   # kinematic viscosity air ~20 degC m2/s


def friction_factor(Re: float, rel_rough: float) -> float:
    """Darcy friction factor: laminar 64/Re, else Colebrook solved iteratively."""
    if Re < 2300:
        return 64.0 / Re
    f = 0.25 / (math.log10(rel_rough / 3.7 + 5.74 / Re ** 0.9)) ** 2  # Swamee-Jain start
    for _ in range(20):
        f = (-2.0 * math.log10(rel_rough / 3.7 + 2.51 / (Re * math.sqrt(f)))) ** -2
    return f


def dp_per_m(flow: float, D: float, rho: float = 1.2, nu: float = AIR_NU,
             eps: float = ROUGHNESS["galv_steel"]) -> float:
    """Pressure gradient [Pa/m] in a round conduit of diameter D [m]."""
    A = math.pi * D ** 2 / 4
    V = flow / A
    f = friction_factor(V * D / nu, eps / D)
    return f / D * rho * V ** 2 / 2


def size_round_for_friction(flow: float, target_pa_m: float, **kw) -> float:
    """Exact diameter [m] giving target friction (bisection, dp decreases with D)."""
    return _bisect(lambda D: target_pa_m - dp_per_m(flow, D, **kw), 0.02, 5.0, tol=1e-7)


def equiv_diameter_rect(a: float, b: float) -> float:
    """Huebscher equivalent round diameter (equal flow & friction)."""
    return 1.30 * (a * b) ** 0.625 / (a + b) ** 0.25


def rect_for_round(D: float, height: float, step: float = 0.05) -> float:
    """Smallest width (multiple of step) at given height with De >= D."""
    w = step
    while equiv_diameter_rect(w, height) < D:
        w += step
    return round(w, 3)


@dataclass
class DuctSection:
    name: str
    flow: float           # m3/s
    length: float         # m
    sum_C: float = 0.0    # sum of fitting loss coefficients (based on section velocity)
    upstream: Optional[str] = None


def size_duct_network(sections: List[DuctSection], target_pa_m: float = 1.0,
                      v_max: float = 8.0, rho: float = 1.2):
    """Equal-friction sizing. Picks standard round size >= exact. Returns rows and
    pressure loss of each path root->leaf (critical path = max)."""
    rows = {}
    for s in sections:
        d_exact = size_round_for_friction(s.flow, target_pa_m, rho=rho)
        # also respect velocity limit
        d_vel = math.sqrt(4 * s.flow / (math.pi * v_max))
        d_need = max(d_exact, d_vel)
        d_std = next((d for d in ROUND_DUCT_MM if d / 1000 >= d_need), ROUND_DUCT_MM[-1]) / 1000
        V = s.flow / (math.pi * d_std ** 2 / 4)
        R = dp_per_m(s.flow, d_std, rho=rho)
        dp = R * s.length + s.sum_C * rho * V ** 2 / 2
        rows[s.name] = dict(sec=s, D=d_std, V=V, R=R, dp=dp)
    # path losses
    children = {s.upstream for s in sections}
    leaves = [s.name for s in sections if s.name not in children]
    paths = {}
    for leaf in leaves:
        tot, n = 0.0, leaf
        while n:
            tot += rows[n]["dp"]
            n = rows[n]["sec"].upstream
        paths[leaf] = tot
    return rows, paths


# ============================================================================
# 8. PIPE SIZING (water)
# ============================================================================
# Steel pipe Sch.40 internal diameter [mm]
STEEL_SCH40_ID = {15: 15.8, 20: 20.9, 25: 26.6, 32: 35.1, 40: 40.9, 50: 52.5, 65: 62.7,
                  80: 77.9, 100: 102.3, 125: 128.2, 150: 154.1, 200: 202.7, 250: 254.5,
                  300: 303.2, 350: 333.4, 400: 381.0}


def water_density(t: float) -> float:
    """kg/m3, Thiesen formula (0..40 degC)."""
    return 1000 * (1 - (t + 288.9414) / (508929.2 * (t + 68.12963)) * (t - 3.9863) ** 2)


def water_viscosity(t: float) -> float:
    """Dynamic viscosity Pa.s (Vogel)."""
    return 2.414e-5 * 10 ** (247.8 / (t + 273.15 - 140))


def water_flow(q_w: float, dT: float, t_mean: float = 9.0) -> float:
    """Volumetric flow [m3/s] for load q [W] and dT [K]. cp = 4190 J/kgK."""
    return q_w / (water_density(t_mean) * 4190 * dT)


def size_pipe(flow: float, t_mean: float = 9.0, max_pa_m: float = 400.0, v_max: float = 3.0,
              eps: float = ROUGHNESS["steel_pipe"]):
    """Smallest Sch.40 DN with friction <= max_pa_m and velocity <= v_max."""
    rho, mu = water_density(t_mean), water_viscosity(t_mean)
    for dn, id_mm in STEEL_SCH40_ID.items():
        D = id_mm / 1000
        V = flow / (math.pi * D ** 2 / 4)
        R = dp_per_m(flow, D, rho=rho, nu=mu / rho, eps=eps)
        if R <= max_pa_m and V <= v_max:
            return dict(DN=dn, ID_mm=id_mm, V=V, R=R)
    raise ValueError("Flow too large for table")


# ============================================================================
# 9. PLANT & POWER
# ============================================================================
KW_PER_TR = 3.51685   # 1 refrigeration ton = 3.517 kW


def fan_power(flow_m3s: float, dp_pa: float, eta_fan=0.70, eta_motor=0.92, eta_drive=0.97):
    """Shaft-to-electric input power [W] = Q dP / (eta_f eta_m eta_d)."""
    return flow_m3s * dp_pa / (eta_fan * eta_motor * eta_drive)


def pump_power(flow_m3s: float, head_pa: float, eta_pump=0.75, eta_motor=0.93):
    return flow_m3s * head_pa / (eta_pump * eta_motor)


def fan_laws(Q1, dp1, P1, n1, n2):
    """Affinity laws: Q ~ n, dP ~ n^2, P ~ n^3."""
    r = n2 / n1
    return Q1 * r, dp1 * r ** 2, P1 * r ** 3


def plant_sizing(coil_total_w: float, cop: float = 5.5, safety: float = 1.10,
                 chw_dT: float = 5.5, cw_dT: float = 5.0, n_chillers: int = 2):
    cap = coil_total_w * safety
    rej = cap * (1 + 1 / cop)
    return dict(chiller_total_kW=cap / 1000, chiller_total_TR=cap / 1000 / KW_PER_TR,
                per_chiller_kW=cap / 1000 / n_chillers, compressor_kW=cap / cop / 1000,
                tower_rejection_kW=rej / 1000,
                chw_flow_Ls=water_flow(cap, chw_dT, 9.0) * 1000,
                cw_flow_Ls=water_flow(rej, cw_dT, 32.0) * 1000)


# ============================================================================
# 10. EQUIPMENT SELECTION & BILL OF MATERIALS (procurement)
# ============================================================================
IEC_MOTOR_KW = [0.18, 0.25, 0.37, 0.55, 0.75, 1.1, 1.5, 2.2, 3, 4, 5.5, 7.5, 11, 15, 18.5,
                22, 30, 37, 45, 55, 75, 90, 110, 132, 160, 200, 250, 315]
KV_SERIES = [0.25, 0.4, 0.63, 1.0, 1.6, 2.5, 4.0, 6.3, 10, 16, 25, 40, 63, 100, 160, 250, 400, 630]
DIFFUSER_NECK_MM = [150, 200, 250, 300, 350, 400]           # round neck, square ceiling diffuser
# Galvanised sheet thickness by largest duct dimension (SMACNA / DW144 style, low pressure)
SHEET_GAUGE = [(300, 0.5), (750, 0.6), (1500, 0.8), (2100, 1.0), (9e9, 1.2)]
# CHW pipe insulation, closed-cell elastomeric, humid climate (condensation control)
CHW_INSUL_MM = [(50, 25), (150, 32), (300, 40), (9e9, 50)]


def motor_size(shaft_kw: float, margin: float = 1.15) -> float:
    """Next standard IEC motor rating >= shaft power x margin."""
    return next(m for m in IEC_MOTOR_KW if m >= shaft_kw * margin)


def control_valve_kv(flow_m3h: float, dp_valve_kpa: float):
    """Kv = Q / sqrt(dp[bar]). Picks nearest standard Kv (log scale), returns
    (Kv_calc, Kv_std, actual dp kPa)."""
    kv = flow_m3h / math.sqrt(dp_valve_kpa / 100)
    std = min(KV_SERIES, key=lambda k: abs(math.log(k / kv)))
    return kv, std, (flow_m3h / std) ** 2 * 100


def select_diffusers(flow_m3s: float, floor_area: float, m2_per_diffuser: float = 16.0,
                     v_neck: float = 3.0):
    """Count by coverage (one per ~4x4 m) and by max ~150 L/s each; neck by velocity."""
    n = max(math.ceil(floor_area / m2_per_diffuser), math.ceil(flow_m3s / 0.150))
    q = flow_m3s / n
    neck = next((d for d in DIFFUSER_NECK_MM if q / (math.pi * (d / 1000) ** 2 / 4) <= v_neck),
                DIFFUSER_NECK_MM[-1])
    return n, q, neck


def select_return_grilles(flow_m3s: float, v_face: float = 2.0, size=(0.6, 0.6)):
    a = size[0] * size[1]
    return math.ceil(flow_m3s / (a * v_face)), size


def duct_sheet_metal(rows: dict, waste: float = 0.15):
    """Surface area, thickness and mass of galvanised round duct sections from
    size_duct_network() rows. Returns list of (name, D_mm, area m2, t mm, kg)."""
    out = []
    for n, r in rows.items():
        d_mm = r["D"] * 1000
        area = math.pi * r["D"] * r["sec"].length * (1 + waste)
        t = next(th for lim, th in SHEET_GAUGE if d_mm <= lim)
        out.append((n, d_mm, area, t, area * t / 1000 * 7850))
    return out


def expansion_tank(system_volume_l: float, t_min: float = 6.0, t_max: float = 35.0,
                   static_bar: float = 3.0, relief_bar: float = 6.0):
    """Diaphragm expansion tank. Ve = V (rho_min/rho_max - 1);
    Vt = Ve / (1 - (Pi+1)/(Pf+1)), Pi = static + 0.3 bar, Pf = 0.9 relief (gauge)."""
    ve = system_volume_l * (water_density(t_min) / water_density(t_max) - 1)
    pi, pf = static_bar + 0.3, 0.9 * relief_bar
    return ve, ve / (1 - (pi + 1) / (pf + 1))


def chw_insulation(dn: int) -> int:
    return next(t for lim, t in CHW_INSUL_MM if dn <= lim)


@dataclass
class BomItem:
    group: str
    item: str
    spec: str
    qty: float
    unit: str
    note: str = ""


def zone_bom(name: str, floor_area: float, V_sa: float, Voz_ls: float, coil_total: float,
             coil_sens: float, on_coil: AirState, off_coil: AirState, condensate_kg_s: float,
             duct_rows: dict, duct_paths: dict, fan_tsp: float, pipe_run_m: float = 15.0,
             coil_dp_kpa: float = 40.0, chw=(6.0, 12.0)) -> List[BomItem]:
    """Bill of materials for one air-handled zone (AHU/ducted FCU + air distribution)."""
    B = []
    dT = chw[1] - chw[0]
    q_w = water_flow(coil_total, dT, sum(chw) / 2)
    esp = max(duct_paths.values()) + 25
    fan_in = fan_power(V_sa, fan_tsp)
    face = V_sa / 2.5
    sel = design_coil(on_coil, V_sa, coil_total, off_coil.t + 0.3, tw_in=chw[0], dT_w=dT)
    coil_txt = (f"{sel[0].rows} rows x {sel[0].fpi:.0f} FPI, 3/8in Cu/Al, {sel[0].circuits} circuits, "
                f"water dp {sel[1]['dp_water_kPa']:.0f} kPa, air dp ~{sel[1]['dp_air']:.0f} Pa"
                if sel else "6-8 rows (rate with design_coil)")
    B.append(BomItem("Air handling", f"AHU / ducted FCU ({name})",
        f"Air {V_sa*1000:.0f} L/s ({V_sa*3600:.0f} m3/h); ESP {esp:.0f} Pa; fan TSP {fan_tsp:.0f} Pa; "
        f"EC/VSD motor {motor_size(fan_in/1000*0.92, 1.25)} kW; CHW coil {coil_txt}, face {face:.2f} m2 @<=2.5 m/s, "
        f"total {coil_total/1000:.1f} kW / sens {coil_sens/1000:.1f} kW; EAT {on_coil.t:.1f}/{on_coil.twb:.1f} WB, "
        f"LAT {off_coil.t:.1f}/{off_coil.twb:.1f} WB; CHW {chw[0]:.0f}/{chw[1]:.0f} C {q_w*1000:.2f} L/s, "
        f"water dp <= {coil_dp_kpa:.0f} kPa; filters G4 + F7; drain pan SS304", 1, "set"))
    n_d, q_d, neck = select_diffusers(V_sa, floor_area)
    B.append(BomItem("Air terminals", "Square ceiling diffuser 4-way, with plenum box & damper",
        f"Neck O{neck} mm, {q_d*1000:.0f} L/s each, NC <= 35, face 600x600", n_d, "pcs"))
    n_r, sz = select_return_grilles(V_sa * 0.9)
    B.append(BomItem("Air terminals", "Return grille (eggcrate / linear) with filter",
        f"{sz[0]*1000:.0f}x{sz[1]*1000:.0f} mm, face velocity <= 2.0 m/s", n_r, "pcs"))
    d_oa = next(d for d in ROUND_DUCT_MM if Voz_ls / 1000 / (math.pi * (d/1000) ** 2 / 4) <= 4.0)
    B.append(BomItem("Air terminals", "Outdoor-air intake louvre + motorised OA damper",
        f"{Voz_ls:.0f} L/s; louvre net free velocity <= 2.5 m/s; OA duct O{d_oa} mm; "
        f"damper actuator 24 VAC 0-10 V", 1, "set"))
    sheets = duct_sheet_metal(duct_rows)
    for t in sorted({s[3] for s in sheets}):
        a = sum(s[2] for s in sheets if s[3] == t)
        kg = sum(s[4] for s in sheets if s[3] == t)
        B.append(BomItem("Ductwork", f"Galvanised steel sheet {t} mm (Z275)",
            f"for round/rect duct, sizes {', '.join(f'O{s[1]:.0f}' for s in sheets if s[3]==t)}; "
            f"incl. 15 % seams/waste", round(kg, 1), "kg", f"{a:.1f} m2"))
    a_tot = sum(s[2] for s in sheets) / 1.15
    B.append(BomItem("Ductwork", "Duct insulation (supply)",
        "Glass wool 25 mm 24-32 kg/m3 with FSK facing, or elastomeric 19 mm", round(a_tot * 1.1, 1), "m2"))
    B.append(BomItem("Ductwork", "Flexible insulated duct to diffusers",
        f"O{neck} mm, <= 1.5 m each, fully extended", n_d * 1.5, "m"))
    n_br = sum(1 for r in duct_rows.values() if r["sec"].upstream)
    B.append(BomItem("Ductwork", "Volume control damper (opposed blade)",
        "one per branch, locking quadrant", n_br, "pcs"))
    B.append(BomItem("Ductwork", "Supports: threaded rod M8/M10, angle 30x30x3, anchors",
        "spacing <= 2.4 m (SMACNA)", math.ceil(sum(r['sec'].length for r in duct_rows.values()) / 2.4),
        "sets"))
    B.append(BomItem("Ductwork", "Flexible canvas connector at AHU", "fire-retardant, 100-150 mm", 2, "pcs"))
    p = size_pipe(q_w, max_pa_m=300)
    ins = chw_insulation(p["DN"])
    B.append(BomItem("Piping", f"CHW branch pipe DN{p['DN']} (steel Sch.40 or PP-R/copper equiv.)",
        f"{q_w*1000:.2f} L/s, V {p['V']:.2f} m/s, R {p['R']:.0f} Pa/m; supply + return",
        2 * pipe_run_m, "m"))
    B.append(BomItem("Piping", "Pipe insulation, closed-cell elastomeric", f"{ins} mm, for DN{p['DN']}",
                     2 * pipe_run_m, "m"))
    trap = math.ceil((fan_tsp * 0.65) / 9.81 + 25)
    B.append(BomItem("Piping", "Condensate drain pipe uPVC DN25-32 + P-trap",
        f"{condensate_kg_s*3600:.1f} kg/h; slope >= 1 %; trap seal >= {trap} mm (negative pressure)",
        pipe_run_m, "m"))
    kv, kvs, dpa = control_valve_kv(q_w * 3600, coil_dp_kpa)
    B.append(BomItem("Valves", "2-way modulating control valve + actuator (or PICV)",
        f"DN{min(p['DN'], 50)}, Kv calc {kv:.2f} -> Kvs {kvs}, dp {dpa:.0f} kPa, equal-%, 24 VAC 0-10 V, "
        f"authority ~{dpa/(dpa+coil_dp_kpa):.2f}", 1, "set"))
    for item, q in [("Isolating ball/butterfly valve", 2), ("Y-strainer 20 mesh", 1),
                    ("Balancing valve with test points (DRV)", 1), ("Flexible connector (rubber)", 2),
                    ("Pressure gauge 0-10 bar + cock", 2), ("Thermometer 0-50 C", 2),
                    ("Automatic air vent", 1), ("Drain valve", 1)]:
        B.append(BomItem("Valves", item, f"DN{p['DN']}, PN16", q, "pcs"))
    for item, spec in [("Zone temperature + RH sensor", "NTC10k / 0-10 V, +-0.3 K, +-3 %RH"),
                       ("CO2 sensor (DCV)", "NDIR 0-2000 ppm, 0-10 V"),
                       ("Supply-air temperature sensor (duct)", "NTC/PT1000, 200 mm probe"),
                       ("Filter differential pressure switch", "50-500 Pa"),
                       ("Fan status (current switch / DP switch)", "DI"),
                       ("Frost/condensate overflow switch", "float switch in drain pan"),
                       ("DDC controller (BACnet MS/TP or IP)", "~8 AI, 6 AO, 6 DI, 4 DO")]:
        B.append(BomItem("Controls", item, spec, 1, "pcs"))
    return B


def plant_bom(coil_total_w: float, n_zones: int, cop: float = 5.5, n_chillers: int = 2,
              t_wb: float = 28.0, chw_head_kpa: float = 250.0, cw_head_kpa: float = 200.0,
              litres_per_kw: float = 8.0, static_bar: float = 3.0) -> List[BomItem]:
    """Central chilled-water plant BOM for n_zones identical zones."""
    P = plant_sizing(coil_total_w * n_zones, cop=cop, n_chillers=n_chillers)
    B = []
    q1 = P["per_chiller_kW"]
    chw1 = P["chw_flow_Ls"] / n_chillers
    cw1 = P["cw_flow_Ls"] / n_chillers
    B.append(BomItem("Plant", "Water-cooled chiller (screw / centrifugal)",
        f"{q1:.0f} kW ({q1/KW_PER_TR:.0f} TR) each; CHW 12->6 C {chw1:.1f} L/s; CW 32->37 C {cw1:.1f} L/s; "
        f"full-load COP >= {cop}, IPLV >= {cop*1.25:.1f}; refrigerant R134a/R513A/R1234ze; "
        f"evap & cond dp <= 70 kPa; 400 V 3ph 50 Hz", n_chillers, "sets"))
    rej1 = P["tower_rejection_kW"] / n_chillers
    B.append(BomItem("Plant", "Cooling tower (induced draft, counter-flow, FRP)",
        f"{rej1:.0f} kW rejection each; {cw1:.1f} L/s; 37->32 C at {t_wb:.0f} C WB "
        f"(approach {32-t_wb:.0f} K); VSD fan; drift <= 0.005 %; SS basin strainer",
        n_chillers, "sets"))
    for nm, fl, head in [("Chilled-water pump", chw1, chw_head_kpa), ("Condenser-water pump", cw1, cw_head_kpa)]:
        shaft = fl / 1000 * head * 1000 / 0.75 / 1000
        B.append(BomItem("Plant", f"{nm} (end-suction / split-case)",
            f"{fl:.1f} L/s ({fl*3.6:.0f} m3/h) @ {head:.0f} kPa ({head/9.81:.0f} m); eta >= 75 %; "
            f"motor {motor_size(shaft)} kW IE3, VSD", n_chillers + 1, "sets", "N duty + 1 standby"))
    for nm, fl in [("CHW", P["chw_flow_Ls"]), ("CW", P["cw_flow_Ls"])]:
        p = size_pipe(fl / 1000, t_mean=9 if nm == "CHW" else 32)
        B.append(BomItem("Plant piping", f"{nm} header pipe DN{p['DN']} steel Sch.40 (welded/grooved)",
            f"{fl:.1f} L/s, V {p['V']:.2f} m/s, R {p['R']:.0f} Pa/m", 0, "m", "length from layout"))
    vol = P["chiller_total_kW"] * litres_per_kw
    ve, vt = expansion_tank(vol, static_bar=static_bar)
    B.append(BomItem("Plant", "Diaphragm expansion tank (CHW)",
        f"system ~{vol:.0f} L; expansion {ve:.0f} L; tank >= {vt:.0f} L, PN10, pre-charge {static_bar+0.3:.1f} bar",
        1, "pcs"))
    for item, spec in [("Air & dirt separator", "on CHW return, header size"),
                       ("Chemical dosing pot + water treatment", "CHW closed loop + CW (inhibitor, biocide, TDS blowdown)"),
                       ("Make-up water tank & transfer pump", f"tower make-up ~{P['tower_rejection_kW']*1.5:.0f} L/h evaporation"),
                       ("Flow meter (electromagnetic) + BTU meter", "CHW main"),
                       ("Chiller plant manager / BMS controller", "sequencing, CHW & CW reset, DP control, BACnet IP"),
                       ("Remote differential pressure sensor (CHW)", "0-250 kPa, at index circuit"),
                       ("Butterfly valves, check valves, strainers, flex connectors (plant)", "per P&ID, PN16")]:
        B.append(BomItem("Plant", item, spec, 1, "lot"))
    B.sort(key=lambda b: b.group == "Plant piping")
    return B


def bom_to_csv(items: List[BomItem], path: str):
    """Write BOM to CSV (opens in Excel). Add your unit prices in the 'unit_price' column."""
    import csv
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["group", "item", "specification", "qty", "unit", "note", "unit_price", "amount"])
        for i, b in enumerate(items, start=2):
            w.writerow([b.group, b.item, b.spec, b.qty, b.unit, b.note, "", f"=D{i}*G{i}"])


def print_bom(items: List[BomItem]):
    g0 = None
    for b in items:
        if b.group != g0:
            print(f"  -- {b.group} --")
            g0 = b.group
        print(f"   {b.qty:>8g} {b.unit:5s} {b.item}\n                  {b.spec}" + (f"  [{b.note}]" if b.note else ""))


# ============================================================================
# 11. HEAT TRANSFER & EQUIPMENT ENGINEERING (component design / manufacturing)
# ============================================================================
AIR_K, AIR_MU, AIR_PR, AIR_CP = 0.0257, 1.81e-5, 0.71, 1006.0   # ~20 degC
WATER_K = 0.58                                                   # W/mK ~10 degC


def lmtd(dt1: float, dt2: float) -> float:
    """Log-mean temperature difference."""
    if abs(dt1 - dt2) < 1e-9:
        return dt1
    return (dt1 - dt2) / math.log(dt1 / dt2)


def eff_counterflow(ntu: float, cr: float) -> float:
    if abs(1 - cr) < 1e-9:
        return ntu / (1 + ntu)
    e = math.exp(-ntu * (1 - cr))
    return (1 - e) / (1 - cr * e)


def eff_crossflow_unmixed(ntu: float, cr: float) -> float:
    """Both fluids unmixed (approximate closed form, Incropera)."""
    if cr < 1e-9:
        return 1 - math.exp(-ntu)
    return 1 - math.exp((1 / cr) * ntu ** 0.22 * (math.exp(-cr * ntu ** 0.78) - 1))


def nu_tube(Re: float, Pr: float, f: Optional[float] = None) -> float:
    """Nusselt in a round tube: laminar 3.66, turbulent Gnielinski (3000<Re<5e6)."""
    if Re < 2300:
        return 3.66
    f = f or (0.790 * math.log(Re) - 1.64) ** -2          # Petukhov smooth tube
    nu = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * math.sqrt(f / 8) * (Pr ** (2 / 3) - 1))
    if Re < 3000:   # linear blend in transition
        return 3.66 + (nu - 3.66) * (Re - 2300) / 700
    return nu


def fin_efficiency_schmidt(h: float, k_fin: float, t_fin: float, r_o: float,
                           Pt: float, Pl: float) -> float:
    """Plate fin on staggered tubes (Schmidt hexagonal equivalent circular fin).
    r_o = collar outer radius [m]; Pt, Pl transverse/longitudinal pitch [m]."""
    XM = Pt / 2
    XL = 0.5 * math.sqrt((Pt / 2) ** 2 + Pl ** 2)
    ratio = 1.27 * (XM / r_o) * math.sqrt(XL / XM - 0.3)
    phi = (ratio - 1) * (1 + 0.35 * math.log(ratio))
    m = math.sqrt(2 * h / (k_fin * t_fin))
    x = m * r_o * phi
    return math.tanh(x) / x


def h_sat_air(t: float, p: float = 101.325) -> float:
    """Enthalpy of saturated air at t [kJ/kg]."""
    return enthalpy(t, w_sat(t, p))


def t_sat_from_h(h: float, p: float = 101.325) -> float:
    return _bisect(lambda t: h_sat_air(t, p) - h, -20.0, 60.0)


@dataclass
class CoilGeometry:
    """Plate-fin round-tube coil. Defaults: 3/8" Cu tube, 25.4 x 22 mm staggered, Al fins."""
    face_w: float = 0.90          # finned length (tube length) m
    face_h: float = 0.381         # height m (= tubes_per_row * Pt)
    rows: int = 6
    fpi: float = 12.0             # fins per inch
    tube_od: float = 9.52e-3
    tube_wall: float = 0.35e-3
    Pt: float = 25.4e-3
    Pl: float = 22.0e-3
    fin_t: float = 0.115e-3
    fin_k: float = 204.0          # Al 1100/8011; Cu = 386
    circuits: int = 15
    fouling_i: float = 0.000088   # m2K/W water side

    @property
    def tubes_per_row(self): return round(self.face_h / self.Pt)


def coil_areas(g: CoilGeometry):
    s_p = 0.0254 / g.fpi
    Dc = g.tube_od + 2 * g.fin_t
    n_fin = g.face_w / s_p
    nt = g.tubes_per_row * g.rows
    depth = g.rows * g.Pl
    A_fin = 2 * (g.face_h * depth - nt * math.pi * Dc ** 2 / 4) * n_fin
    A_tube = math.pi * Dc * (g.face_w - n_fin * g.fin_t) * nt
    A_o = A_fin + A_tube
    A_fr = g.face_w * g.face_h
    A_min = A_fr * ((g.Pt - Dc) / g.Pt) * ((s_p - g.fin_t) / s_p)
    Di = g.tube_od - 2 * g.tube_wall
    A_i = math.pi * Di * g.face_w * nt
    return dict(s_p=s_p, Dc=Dc, A_fin=A_fin, A_o=A_o, A_fr=A_fr, A_min=A_min, Di=Di, A_i=A_i,
                n_fin=n_fin, nt=nt, depth=depth)


def h_air_gray_webb(g: CoilGeometry, A: dict, m_air: float):
    """Air-side h [W/m2K] - Gray & Webb (1986) plain plate-fin correlation."""
    Gc = m_air / A["A_min"]
    Re = Gc * A["Dc"] / AIR_MU
    s = A["s_p"] - g.fin_t
    j4 = 0.14 * Re ** -0.328 * (g.Pt / g.Pl) ** -0.502 * (s / A["Dc"]) ** 0.0312
    j = j4
    if g.rows < 4:
        j = j4 * 0.991 * (2.24 * Re ** -0.092 * (g.rows / 4) ** -0.031) ** (0.607 * (4 - g.rows))
    return j * Gc * AIR_CP / AIR_PR ** (2 / 3), Re


def rate_coil(g: CoilGeometry, air_in: AirState, V_air: float, tw_in: float, q_water: float):
    """Rate a chilled-water cooling coil (counter-flow circuiting).
    Wet surface: Braun-Klein-Mitchell enthalpy effectiveness model.
    Dry surface: classic epsilon-NTU. Returns dict with capacity and leaving states."""
    p = air_in.p
    A = coil_areas(g)
    m_a = V_air / air_in.v                         # kg dry air/s
    ho, Re_a = h_air_gray_webb(g, A, m_a)
    # water side
    tw_m = tw_in + 3
    rho_w, mu_w = water_density(tw_m), water_viscosity(tw_m)
    m_w = q_water * rho_w
    v_w = q_water / (g.circuits * math.pi * A["Di"] ** 2 / 4)
    Re_w = rho_w * v_w * A["Di"] / mu_w
    Pr_w = mu_w * 4190 / WATER_K
    hi = nu_tube(Re_w, Pr_w) * WATER_K / A["Di"]
    r_o = A["Dc"] / 2
    face_v = V_air / A["A_fr"]

    # --- try wet-coil model, iterate on cs (saturation-enthalpy slope) ---
    tw_out = tw_in + 5.0
    for _ in range(30):
        cs = (h_sat_air(tw_out, p) - h_sat_air(tw_in, p)) / max(tw_out - tw_in, 0.05) * 1000  # J/kgK
        eta_f = fin_efficiency_schmidt(ho * cs / AIR_CP, g.fin_k, g.fin_t, r_o, g.Pt, g.Pl)
        eta_o = 1 - A["A_fin"] / A["A_o"] * (1 - eta_f)
        R = cs / (hi * A["A_i"]) + cs * g.fouling_i / A["A_i"] + AIR_CP / (ho * eta_o * A["A_o"])
        C_a, C_w = m_a, m_w * 4190 / cs
        Cmin, Cmax = min(C_a, C_w), max(C_a, C_w)
        ntu = 1 / (R * Cmin)
        eps = eff_counterflow(ntu, Cmin / Cmax)
        Q = eps * Cmin * (air_in.h - h_sat_air(tw_in, p)) * 1000                  # W
        new = tw_in + Q / (m_w * 4190)
        if abs(new - tw_out) < 1e-4:
            break
        tw_out = 0.5 * (tw_out + new)
    h_out = air_in.h - Q / m_a / 1000
    ntu_o = ho * eta_o * A["A_o"] / (m_a * AIR_CP)
    hs_eff = air_in.h - (air_in.h - h_out) / (1 - math.exp(-ntu_o))
    ts_eff = t_sat_from_h(hs_eff, p)
    t_out = ts_eff + (air_in.t - ts_eff) * math.exp(-ntu_o)
    w_out = (h_out - CP_DA * t_out) / (HFG0 + CP_V * t_out)
    wet = ts_eff < air_in.tdp and w_out < air_in.w
    if not wet:   # dry coil: sensible epsilon-NTU
        eta_f = fin_efficiency_schmidt(ho, g.fin_k, g.fin_t, r_o, g.Pt, g.Pl)
        eta_o = 1 - A["A_fin"] / A["A_o"] * (1 - eta_f)
        UA = 1 / (1 / (hi * A["A_i"]) + g.fouling_i / A["A_i"] + 1 / (ho * eta_o * A["A_o"]))
        Ca, Cw = m_a * (CP_DA + CP_V * air_in.w) * 1000, m_w * 4190
        Cmin, Cmax = min(Ca, Cw), max(Ca, Cw)
        eps = eff_counterflow(UA / Cmin, Cmin / Cmax)
        Q = eps * Cmin * (air_in.t - tw_in)
        t_out, w_out = air_in.t - Q / Ca, air_in.w
        tw_out = tw_in + Q / Cw
        h_out = enthalpy(t_out, w_out)
    out = AirState(t_out, w_out, p)
    Q_s = m_a * (CP_DA + CP_V * air_in.w) * (air_in.t - t_out) * 1000
    # water pressure drop: straight tubes + return bends (K=1.5 each)
    tubes_per_circ = A["nt"] / g.circuits
    L = g.face_w * tubes_per_circ
    f = friction_factor(max(Re_w, 1), 1.5e-6 / A["Di"])
    dpw = (f * L / A["Di"] + 1.5 * (tubes_per_circ - 1)) * rho_w * v_w ** 2 / 2
    # air pressure drop: EMPIRICAL estimate (+-30 %), wet x1.4 - verify with maker's software
    dpa = g.rows * 3.6 * (g.fpi / 12) ** 0.5 * face_v ** 1.75 * (1.4 if wet else 1.0)
    return dict(Q=Q, Q_s=Q_s, SHR=Q_s / Q if Q else 1, air_out=out, tw_out=tw_out, wet=wet,
                face_v=face_v, v_water=v_w, Re_air=Re_a, Re_water=Re_w, h_o=ho, h_i=hi,
                eta_fin=eta_f, A_o=A["A_o"], A_i=A["A_i"], dp_water_kPa=dpw / 1000, dp_air=dpa,
                condensate_kg_h=m_a * max(0, air_in.w - w_out) * 3600, ntu=ntu if wet else UA / Cmin)


def design_coil(air_in: AirState, V_air: float, Q_req: float, t_leave_max: float,
                tw_in: float = 6.0, dT_w: float = 6.0, face_v_max: float = 2.5,
                rows_opts=(4, 6, 8), fpi_opts=(10, 12, 14)):
    """Smallest coil (rows x fpi) meeting capacity and leaving-air DB at face velocity
    <= face_v_max. Face is sized from tube pitch; water circuits for 0.6-1.5 m/s."""
    q_w = Q_req / (water_density(tw_in + dT_w / 2) * 4190 * dT_w)
    best = None
    for rows in rows_opts:
        for fpi in fpi_opts:
            g = CoilGeometry(rows=rows, fpi=fpi)
            nt = math.ceil(V_air / face_v_max / g.face_w / g.Pt)
            g.face_h = nt * g.Pt
            Di = g.tube_od - 2 * g.tube_wall
            g.circuits = max(1, min(nt, round(q_w / (1.0 * math.pi * Di ** 2 / 4))))
            r = rate_coil(g, air_in, V_air, tw_in, q_w)
            ok = r["Q"] >= Q_req and r["air_out"].t <= t_leave_max
            if ok and (best is None or rows * fpi < best[0].rows * best[0].fpi):
                best = (g, r)
    return best


def insulation_condensation(t_fluid: float, t_amb: float, rh_amb: float, d_out_pipe: float,
                            k_ins: float = 0.036, h_surf: float = 9.0, margin: float = 0.5):
    """Minimum insulation thickness [m] on a cold pipe so the jacket surface stays above
    ambient dew point + margin.  Cylinder: R_ins = ln(r2/r1)/(2 pi k); R_s = 1/(h 2 pi r2)."""
    tdp = dew_point(w_from_rh(t_amb, rh_amb))
    r1 = d_out_pipe / 2

    def surface_t(th):
        r2 = r1 + th
        Ri = math.log(r2 / r1) / (2 * math.pi * k_ins)
        Rs = 1 / (h_surf * 2 * math.pi * r2)
        q = (t_amb - t_fluid) / (Ri + Rs)
        return t_amb - q * Rs
    th = _bisect(lambda th: surface_t(th) - (tdp + margin), 1e-4, 0.3)
    return th, tdp


def npsh_available(p_surface_kpa: float, t_water: float, z_static_m: float,
                   h_friction_m: float) -> float:
    """NPSHa [m] = (p_s - p_v)/(rho g) + z_s - h_f  (z_s + if liquid level above pump)."""
    rho = water_density(t_water)
    return (p_surface_kpa - p_ws(t_water)) * 1000 / (rho * 9.81) + z_static_m - h_friction_m


def vapour_compression(h1: float, h2s: float, h3: float, Q_evap_kw: float,
                       eta_is: float = 0.70, v1: Optional[float] = None, eta_vol: float = 0.90):
    """Single-stage cycle from refrigerant-table enthalpies [kJ/kg]:
    1 = evaporator outlet (superheated), 2s = isentropic discharge, 3 = condenser outlet
    (subcooled), 4 = h3 (isenthalpic expansion). Returns mass flow, power, COP, displacement."""
    h2 = h1 + (h2s - h1) / eta_is
    q_e = h1 - h3
    m = Q_evap_kw / q_e
    W = m * (h2 - h1)
    out = dict(h2=h2, refrigerating_effect=q_e, mass_flow=m, compressor_kW=W,
               condenser_kW=m * (h2 - h3), COP=Q_evap_kw / W)
    if v1:
        out["displacement_m3h"] = m * v1 / eta_vol * 3600
    return out


def db_sum(levels: List[float]) -> float:
    """Logarithmic addition of sound levels: L = 10 log10(sum 10^(Li/10))."""
    return 10 * math.log10(sum(10 ** (L / 10) for L in levels))


def room_spl(Lw: float, r: float, room_const_R: float, Q_dir: float = 2.0) -> float:
    """Sound pressure level in a room: Lp = Lw + 10 log(Q/(4 pi r^2) + 4/R)."""
    return Lw + 10 * math.log10(Q_dir / (4 * math.pi * r ** 2) + 4 / room_const_R)


def vibration_isolation(rpm: float, static_deflection_mm: float):
    """Spring isolator: fn = 15.76/sqrt(delta_mm) Hz; transmissibility T = 1/|r^2 - 1|;
    isolation efficiency = 1 - T (undamped)."""
    fn = 15.76 / math.sqrt(static_deflection_mm)
    r = (rpm / 60) / fn
    T = 1 / abs(r ** 2 - 1)
    return fn, T, (1 - T) if r > math.sqrt(2) else 0.0


def duct_leakage_en12237(area_m2: float, p_pa: float, tclass: str = "B") -> float:
    """Max allowed leakage [L/s] for duct surface area at test pressure.
    f_max = K * p^0.65 [L/(s.m2)], K: A 0.027, B 0.009, C 0.003, D 0.001."""
    K = {"A": 0.027, "B": 0.009, "C": 0.003, "D": 0.001}[tclass]
    return K * p_pa ** 0.65 * area_m2


# ============================================================================
# WORKED EXAMPLE
# ============================================================================
def _example():
    line = "=" * 78
    print(line + "\nHVAC ENGINE - WORKED EXAMPLE : open-plan office, 21.0 N 105.8 E, July 21\n" + line)

    site = Site(latitude=21.0, longitude=105.8, tz=7, altitude=10, tau_b=0.45, tau_d=2.2)
    p = pressure_at_altitude(site.altitude)
    T_DB, DR, T_WB = 36.0, 7.0, 28.0          # outdoor design: DB, daily range, MCWB
    day = 202                                  # July 21

    print("\n[1] Psychrometrics")
    oa = AirState.from_wb(T_DB, T_WB, p)
    ra = AirState.from_rh(24.0, 0.50, p)
    print(f"  Outdoor : {oa}\n  Room    : {ra}")

    print("\n[2] Zone cooling load (24 h)")
    z = Zone("Office-01", floor_area=150, height=3.0, t_in=24, rh_in=0.50, people=15,
             lighting_wpm2=9, equipment_wpm2=12, infiltration_ach=0.3, construction="medium")
    z.walls = [Opaque("W-wall", 30 - 12, 0.6, 90, 90, 0.6, 0.45, 7),
               Opaque("S-wall", 45 - 18, 0.6, 90, 0, 0.6, 0.45, 7),
               Opaque("Roof", 150, 0.45, 0, 0, 0.7, 0.35, 6)]
    z.windows = [Window("W-glass", 12, 2.8, 0.40, 90, shade_coeff=0.8),
                 Window("S-glass", 18, 2.8, 0.40, 0, shade_coeff=0.8)]
    res = cooling_load(z, site, day, T_DB, DR, T_WB)
    print("  hour  T_out  walls  roof  g_cnd  solar  peopl  light  equip  infil | SENS   LAT  TOTAL")
    for h in range(24):
        c = res["components"]
        s, l = res["sensible"][h], res["latent"][h]
        print(f"  {h+1:4d} {res['outdoor'][h].t:6.1f}" +
              "".join(f"{c[k][h]:7.0f}" for k in ["walls", "roof", "glass_cond", "solar", "people",
                                                    "lighting", "equipment", "infiltration_s"]) +
              f" |{s:6.0f}{l:6.0f}{s+l:7.0f}")
    ph = res["peak_hour"]
    Qs, Ql = res["sensible"][ph], res["latent"][ph]
    print(f"  PEAK at hour {ph+1}: sensible {Qs/1000:.2f} kW, latent {Ql/1000:.2f} kW, "
          f"total {(Qs+Ql)/1000:.2f} kW, RSHR {Qs/(Qs+Ql):.2f}, {(Qs+Ql)/z.floor_area:.0f} W/m2")

    print("\n[3] Ventilation (ASHRAE 62.1)")
    Vbz, Voz = outdoor_air("office", z.floor_area, z.people, Ez=1.0)
    print(f"  Vbz = 2.5*{z.people} + 0.3*{z.floor_area} = {Vbz:.1f} L/s ; Voz = {Voz:.1f} L/s")

    print("\n[4] Supply air & coil")
    t_sup = 13.0
    V_sa = supply_airflow(Qs, ra.t, t_sup, rho=ra.rho)
    w_sa = supply_humidity(ra, Ql, V_sa, rho=ra.rho)
    sa = AirState(t_sup, w_sa, p)
    oa_frac = (Voz / 1000) / V_sa
    ma = mix(oa, oa_frac, ra, 1 - oa_frac)
    fan_dT = 0.8                                # draw-through fan heat
    off_coil = AirState(t_sup - fan_dT, w_sa, p)
    adp, bf = coil_adp_bf(ma, off_coil)
    Qc, Qcs, cond = coil_load(ma, off_coil, V_sa)
    print(f"  Supply airflow  : {V_sa:.3f} m3/s = {V_sa*1000:.0f} L/s ({V_sa/z.floor_area*1000:.1f} L/s.m2,"
          f" {V_sa*3600/(z.floor_area*z.height):.1f} ACH)")
    print(f"  OA fraction     : {oa_frac*100:.1f} %")
    print(f"  Supply  : {sa}\n  Mixed   : {ma}\n  Off-coil: {off_coil}")
    print(f"  ADP {adp:.2f} C, bypass factor {bf:.3f}")
    print(f"  Coil load: total {Qc/1000:.2f} kW, sensible {Qcs/1000:.2f} kW, "
          f"condensate {cond*3600:.2f} kg/h, SHR {Qcs/Qc:.2f}")

    print("\n[5] Duct sizing (equal friction 1.0 Pa/m, galv. steel)")
    q = V_sa
    secs = [DuctSection("AHU-A", q, 8, 1.2),
            DuctSection("A-B", q * 0.5, 6, 0.3, "AHU-A"),
            DuctSection("B-D1", q * 0.25, 4, 1.5, "A-B"),
            DuctSection("B-D2", q * 0.25, 7, 1.8, "A-B"),
            DuctSection("A-C", q * 0.5, 10, 0.6, "AHU-A"),
            DuctSection("C-D3", q * 0.25, 3, 1.5, "A-C"),
            DuctSection("C-D4", q * 0.25, 6, 1.8, "A-C")]
    rows, paths = size_duct_network(secs, 1.0, 8.0)
    print("  section   flow L/s   D mm   rect @250mm   V m/s   R Pa/m   dP Pa")
    for n, r in rows.items():
        print(f"  {n:8s} {r['sec'].flow*1000:8.0f} {r['D']*1000:6.0f}"
              f"   {rect_for_round(r['D'], 0.25)*1000:5.0f}x250  {r['V']:7.2f} {r['R']:7.2f} {r['dp']:7.1f}")
    crit = max(paths, key=paths.get)
    print(f"  Critical path -> {crit}: {paths[crit]:.1f} Pa (+ diffuser 25 Pa, coil 150 Pa, filter 120 Pa)")
    esp = paths[crit] + 25
    tsp = esp + 150 + 120 + 50
    print(f"  Fan total pressure ~ {tsp:.0f} Pa -> fan input power {fan_power(V_sa, tsp):.0f} W")

    print("\n[6] Heating load (winter 8 C outside, 22 C inside)")
    elems = [(0.6, 45), (0.45, 150), (2.8, 30)]
    print(f"  Q_heat = {heating_load(elems, 22, 8, 450, 0.5):.0f} W")

    print("\n[7] Plant - building of 40 such zones")
    plant = plant_sizing(Qc * 40, cop=5.5, n_chillers=2)
    for k, v in plant.items():
        print(f"  {k:20s}: {v:9.1f}")
    pipe = size_pipe(plant["chw_flow_Ls"] / 1000)
    print(f"  CHW header: DN{pipe['DN']}  V {pipe['V']:.2f} m/s  R {pipe['R']:.0f} Pa/m")
    head = 250e3
    print(f"  CHW pump (250 kPa): {pump_power(plant['chw_flow_Ls']/1000, head)/1000:.1f} kW")

    print("\n[8] Equipment schedule & bill of materials")
    print("  ZONE (one of 40):")
    zb = zone_bom(z.name, z.floor_area, V_sa, Voz, Qc, Qcs, ma, off_coil, cond, rows, paths, tsp)
    print_bom(zb)
    print("  CENTRAL PLANT (40 zones):")
    pb = plant_bom(Qc, 40, cop=5.5, n_chillers=2, t_wb=T_WB)
    print_bom(pb)
    import os
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bom_example.csv")
    bom_to_csv(zb + pb, out)
    print(f"  CSV written: {out}")

    print("\n[9] Equipment engineering")
    print("  9a Cooling-coil design (rate candidate geometries, pick smallest that meets duty)")
    best = design_coil(ma, V_sa, Qc, off_coil.t + 0.3, tw_in=6.0, dT_w=6.0)
    if best:
        g, r = best
        print(f"   Geometry : {g.face_w*1000:.0f} x {g.face_h*1000:.0f} mm face, {g.rows} rows, {g.fpi:.0f} FPI, "
              f"3/8in Cu tube {g.tubes_per_row} tubes/row, Al fin {g.fin_t*1000:.3f} mm, {g.circuits} circuits")
        print(f"   Areas    : outside {r['A_o']:.1f} m2 (ratio Ao/Ai {r['A_o']/r['A_i']:.1f}), "
              f"face velocity {r['face_v']:.2f} m/s")
        print(f"   Air side : Re {r['Re_air']:.0f}, h_o {r['h_o']:.1f} W/m2K, fin eff {r['eta_fin']:.2f}; "
              f"water: v {r['v_water']:.2f} m/s, Re {r['Re_water']:.0f}, h_i {r['h_i']:.0f} W/m2K")
        print(f"   Rating   : Q {r['Q']/1000:.2f} kW (req {Qc/1000:.2f}), SHR {r['SHR']:.2f}, "
              f"wet={r['wet']}, LAT {r['air_out'].t:.1f} C / {r['air_out'].w*1000:.2f} g/kg, "
              f"water out {r['tw_out']:.1f} C")
        print(f"   dP       : water {r['dp_water_kPa']:.1f} kPa, air ~{r['dp_air']:.0f} Pa (empirical), "
              f"condensate {r['condensate_kg_h']:.1f} kg/h")
    print("  9b Insulation for condensation control (ambient 32 C / 80 %RH, CHW 6 C)")
    for dn in (25, 50, 100, 150, 200):
        od = {25: 33.4, 50: 60.3, 100: 114.3, 150: 168.3, 200: 219.1}[dn] / 1000
        th, tdp = insulation_condensation(6.0, 32.0, 0.80, od, k_ins=0.036)
        print(f"   DN{dn:<4d} OD {od*1000:5.1f} mm -> min {th*1000:4.1f} mm (dew point {tdp:.1f} C)")
    print("  9c Pump NPSH (open tower basin 1.5 m above CW pump, 32 C, suction loss 0.8 m)")
    print(f"   NPSHa = {npsh_available(101.3, 32, 1.5, 0.8):.1f} m  (need NPSHr + 0.5..1.5 m margin)")
    print("  9d Vapour-compression cycle, R134a, Te 0 C (SH 5 K), Tc 40 C (SC 5 K), 100 kW")
    cyc = vapour_compression(h1=402.7, h2s=427.0, h3=249.7, Q_evap_kw=100, eta_is=0.72,
                             v1=0.0709, eta_vol=0.90)
    print("   " + ", ".join(f"{k} {v:.3f}" for k, v in cyc.items()))
    print("  9e Acoustics: AHU Lw 75 dB + 2 x 70 dB sources ->"
          f" {db_sum([75, 70, 70]):.1f} dB; Lp at 3 m, R=200 m2: {room_spl(db_sum([75,70,70]), 3, 200):.1f} dB")
    fn, T, eff = vibration_isolation(1450, 25)
    print(f"  9f Vibration: fan 1450 rpm on 25 mm springs -> fn {fn:.2f} Hz, T {T:.3f}, isolation {eff*100:.1f} %")
    a_duct = sum(s[2] for s in duct_sheet_metal(rows)) / 1.15
    print(f"  9g Duct leakage EN 12237 class B @ 400 Pa, {a_duct:.1f} m2: max "
          f"{duct_leakage_en12237(a_duct, 400, 'B'):.2f} L/s ({duct_leakage_en12237(a_duct,400,'B')/(V_sa*1000)*100:.2f} % of flow)")


if __name__ == "__main__":
    _example()
