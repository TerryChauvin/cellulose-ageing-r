"""Generate data/exemple_vieillissement.csv: SYNTHETIC data, not measurements.

Model: Ekenstam kinetics 1/DP(t) = 1/DP0 + k(T) t, Arrhenius k(T), inverse of
DP^0.905 = 0.75 [eta], multiplicative noise. Fixed seed, stdlib only.
"""
import csv
import math
import random

random.seed(20404)
R = 8.314
# material: (DP0, k at 80 C [1/day], Ea [J/mol], pH0, pH loss/28 d, RA0, RA loss/28 d)
MATERIALS = {
    "carton_tamponne": (1200, 1.5e-5, 105e3, 8.6, 0.3, 0.30, 0.04),
    "carton_acide": (1000, 6.0e-5, 95e3, 5.2, 0.6, 0.00, 0.00),
    "papier_non_tamponne": (900, 3.5e-5, 100e3, 6.8, 0.8, 0.02, 0.01),
}
TEMPS = [70, 80, 90]
DAYS = [0, 3, 7, 14, 28]
RH = 50
REPS = 3

rows = []
for mat, (dp0, k80, ea, ph0, dph, ra0, dra) in MATERIALS.items():
    for temp in TEMPS:
        k = k80 * math.exp(-ea / R * (1 / (temp + 273.15) - 1 / 353.15))
        for rep in range(REPS):
            dp0_i = dp0 * random.gauss(1, 0.03)
            for d in DAYS:
                dp = 1 / (1 / dp0_i + k * d) * random.gauss(1, 0.02)
                eta = dp**0.905 / 0.75
                frac = d / 28 * (k / k80)
                ph = ph0 - dph * frac + random.gauss(0, 0.05)
                ra = max(0.0, ra0 - dra * frac + random.gauss(0, 0.01)) if ra0 > 0 else 0.0
                rows.append([mat, temp, RH, d, round(eta, 1), round(ph, 2), round(ra, 3)])

with open("data/exemple_vieillissement.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["materiau", "temperature_C", "humidite_pct", "jours",
                "viscosite_intrinseque_mL_g", "pH", "reserve_alcaline_mol_kg"])
    w.writerows(rows)
print(len(rows), "rows")
