#!/usr/bin/env python3
"""
nirps_mag_estimator.py
======================
Estimate any magnitude of a main-sequence star using the
Pecaut & Mamajek (2013) dwarf color/Teff table.

Given:
  - Effective temperature (Teff) OR spectral type
  - At least one known magnitude (e.g. J, H, V, G, Ks...)
  - Distance (optional, for absolute magnitude conversion)

Outputs all magnitudes available in the Pecaut-Mamajek table
by interpolating the color indices from the nearest Teff match.

Usage:
    python nirps_mag_estimator.py

Edit the CONFIGURATION section below.

Requirements:
    pip install numpy requests
"""

import numpy as np
import urllib.request
import sys

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit this section
# ══════════════════════════════════════════════════════════════════════════════

TARGET_NAME  = "TOI-4349"

# Provide Teff OR spectral type (Teff takes priority if both given)
TEFF         = 3850              # effective temperature in K (set to None to use SpType)
SPTYPE       = None              # e.g. "F5V", "K2V", "M3V" (used if TEFF is None)

# Known magnitudes — fill in what you have, leave others as None
# At least ONE must be provided
KNOWN_MAGS = {
    "V":   None,    # Johnson V
    "J":   10.4,   # 2MASS J
    "H":   None,   # 2MASS H
    "Ks":  None,    # 2MASS Ks
    "G":   None,    # Gaia G
    "B":   None,    # Johnson B
    "I":   None,    # Cousins Ic
    "W1":  None,    # WISE W1
    "W2":  None,    # WISE W2
}

# Distance (optional) — provide ONE of these for absolute magnitude conversion
PARALLAX_MAS    = 13.89          # parallax in milliarcseconds (e.g. from Gaia)
DISTANCE_PC     = None           # distance in parsecs

# Table URL (fetched online — requires internet connection)
TABLE_URL = "https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt"

# ══════════════════════════════════════════════════════════════════════════════
# Parse Pecaut-Mamajek table
# ══════════════════════════════════════════════════════════════════════════════

# Column definitions matching the table header
# #SpT Teff logT BCv logL Mbol R_Rsun Mv B-V Bt-Vt G-V Bp-Rp G-Rp M_G b-y U-B V-Rc V-Ic V-Ks J-H H-Ks M_J M_Ks Ks-W1 W1-W2 W1-W3 W1-W4 g-r i-z z-Y Msun
COL_NAMES = [
    "SpT","Teff","logT","BCv","logL","Mbol","R_Rsun","Mv",
    "B-V","Bt-Vt","G-V","Bp-Rp","G-Rp","M_G","b-y","U-B",
    "V-Rc","V-Ic","V-Ks","J-H","H-Ks","M_J","M_Ks",
    "Ks-W1","W1-W2","W1-W3","W1-W4","g-r","i-z","z-Y","Msun"
]

def fetch_table(url):
    """
    Download and parse the Pecaut-Mamajek table.
    First tries a local file 'EEM_dwarf_UBVIJHK_colors_Teff.txt',
    then falls back to fetching from the web.
    """
    import os
    local_file = "EEM_dwarf_UBVIJHK_colors_Teff.txt"
    if os.path.exists(local_file):
        print(f"  Reading Pecaut-Mamajek table from local file: {local_file}")
        with open(local_file, "r") as f:
            lines = f.read().splitlines()
    else:
        print(f"  Fetching Pecaut-Mamajek table from web...")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=15) as r:
                lines = r.read().decode("utf-8").splitlines()
        except Exception as e:
            print(f"  ERROR: Could not fetch table: {e}")
            print(f"  Download the table manually from:")
            print(f"  {url}")
            print(f"  Save it as 'EEM_dwarf_UBVIJHK_colors_Teff.txt' in the same folder.")
            sys.exit(1)

    rows = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Only parse lines that start with a spectral type (letter)
        if not line[0].isalpha():
            continue
        # Stop at notes section
        if line.startswith("The adopted") or line.startswith("Further"):
            break
        parts = line.split()
        if len(parts) < 10:
            continue
        row = {"SpT": parts[0]}
        try:
            row["Teff"] = float(parts[1])
        except ValueError:
            continue
        # Parse remaining columns
        for i, col in enumerate(COL_NAMES[2:], start=2):
            try:
                val = parts[i]
                row[col] = float(val) if val not in ("...", ".....", "....", "..") else np.nan
            except (IndexError, ValueError):
                row[col] = np.nan
        rows.append(row)

    print(f"  Parsed {len(rows)} spectral type entries.")
    return rows

def find_nearest_rows(rows, teff):
    """Find the two nearest rows by Teff for interpolation."""
    teffs = np.array([r["Teff"] for r in rows])
    idx   = np.argsort(np.abs(teffs - teff))
    # Return up to 2 nearest rows for interpolation
    i1, i2 = idx[0], idx[1]
    return rows[i1], rows[i2]

def interpolate_row(r1, r2, teff):
    """Linearly interpolate all columns between two rows at given Teff."""
    t1, t2 = r1["Teff"], r2["Teff"]
    if abs(t2 - t1) < 1:
        return r1
    alpha = (teff - t1) / (t2 - t1)
    result = {"Teff": teff}
    for col in COL_NAMES[1:]:
        v1 = r1.get(col, np.nan)
        v2 = r2.get(col, np.nan)
        if np.isnan(v1) or np.isnan(v2):
            result[col] = v1 if not np.isnan(v1) else v2
        else:
            result[col] = v1 + alpha * (v2 - v1)
    return result

def sptype_to_teff(rows, sptype):
    """Look up Teff for a given spectral type string."""
    sptype = sptype.strip().upper()
    for r in rows:
        if r["SpT"].upper() == sptype:
            return r["Teff"]
    # Try partial match
    for r in rows:
        if r["SpT"].upper().startswith(sptype[:3]):
            return r["Teff"]
    return None

def distance_modulus(parallax_mas=None, distance_pc=None):
    """Compute distance modulus from parallax or distance."""
    if parallax_mas is not None and parallax_mas > 0:
        d_pc = 1000.0 / parallax_mas
    elif distance_pc is not None and distance_pc > 0:
        d_pc = distance_pc
    else:
        return None, None
    mu = 5 * np.log10(d_pc) - 5
    return mu, d_pc

# ══════════════════════════════════════════════════════════════════════════════
# Magnitude estimator
# ══════════════════════════════════════════════════════════════════════════════

# Map between band names and the color indices in the PM table
# Each entry: (band, ref_band, color_col) means band = ref_band - color
# Or (band, ref_band, color_col, sign=-1) means band = ref_band + color
BAND_RELATIONS = {
    # From V
    "B":   ("V", "B-V",   +1),   # B  = V + (B-V)
    "U":   ("B", "U-B",   +1),   # U  = B + (U-B)
    "Rc":  ("V", "V-Rc",  -1),   # Rc = V - (V-Rc)
    "Ic":  ("V", "V-Ic",  -1),   # Ic = V - (V-Ic)
    "Ks":  ("V", "V-Ks",  -1),   # Ks = V - (V-Ks)
    "G":   ("V", "G-V",   +1),   # G  = V + (G-V)   [Gaia]
    # From Ks
    "J":   ("Ks", "J-H",  None, "H", "H-Ks"),  # J = Ks + (H-Ks) + (J-H)
    "H":   ("Ks", "H-Ks", +1),   # H  = Ks + (H-Ks)
    "W1":  ("Ks", "Ks-W1",-1),   # W1 = Ks - (Ks-W1)
    "W2":  ("W1", "W1-W2",-1),   # W2 = W1 - (W1-W2)
    "W3":  ("W1", "W1-W3",-1),   # W3 = W1 - (W1-W3)
    "W4":  ("W1", "W1-W4",-1),   # W4 = W1 - (W1-W4)
}

def estimate_all_mags(known_mags, pm_row):
    """
    Given a set of known magnitudes and a PM table row (with color indices),
    estimate all other magnitudes by applying color indices.
    Uses multiple paths and takes the average when multiple estimates exist.
    """
    mags = {k: v for k, v in known_mags.items() if v is not None}

    # First derive V if not known, from any available band
    if "V" not in mags:
        # Try V from J via V-Ks and then Ks-J
        for band, color_col in [("J","V-Ks"), ("H","V-Ks"), ("Ks","V-Ks")]:
            if band in mags and not np.isnan(pm_row.get(color_col, np.nan)):
                # V-Ks known, and band offset to Ks known
                if band == "Ks":
                    mags["V"] = mags["Ks"] + pm_row["V-Ks"]
                elif band == "H":
                    ks_est = mags["H"] - pm_row.get("H-Ks", np.nan)
                    if not np.isnan(ks_est):
                        mags["V"] = ks_est + pm_row["V-Ks"]
                elif band == "J":
                    hks = pm_row.get("H-Ks", np.nan)
                    jh  = pm_row.get("J-H", np.nan)
                    if not np.isnan(hks) and not np.isnan(jh):
                        ks_est = mags["J"] - jh - hks
                        mags["V"] = ks_est + pm_row["V-Ks"]
                if "V" in mags:
                    break

    # Now derive everything from V and Ks using simple color relations
    max_iter = 5
    for _ in range(max_iter):
        prev_len = len(mags)

        # V-based
        if "V" in mags:
            v = mags["V"]
            for col, color_col, sign in [
                ("B",  "B-V",  +1),
                ("Rc", "V-Rc", -1),
                ("Ic", "V-Ic", -1),
                ("Ks", "V-Ks", -1),
                ("G",  "G-V",  +1),
            ]:
                if col not in mags:
                    c = pm_row.get(color_col, np.nan)
                    if not np.isnan(c):
                        mags[col] = v + sign * c

        # U from B
        if "B" in mags and "U" not in mags:
            c = pm_row.get("U-B", np.nan)
            if not np.isnan(c):
                mags["U"] = mags["B"] + c

        # Ks-based
        if "Ks" in mags:
            ks = mags["Ks"]
            if "H" not in mags:
                c = pm_row.get("H-Ks", np.nan)
                if not np.isnan(c):
                    mags["H"] = ks + c
            if "J" not in mags and "H" in mags:
                c = pm_row.get("J-H", np.nan)
                if not np.isnan(c):
                    mags["J"] = mags["H"] + c
            if "W1" not in mags:
                c = pm_row.get("Ks-W1", np.nan)
                if not np.isnan(c):
                    mags["W1"] = ks - c

        # W1-based
        if "W1" in mags:
            w1 = mags["W1"]
            for col, color_col in [("W2","W1-W2"),("W3","W1-W3"),("W4","W1-W4")]:
                if col not in mags:
                    c = pm_row.get(color_col, np.nan)
                    if not np.isnan(c):
                        mags[col] = w1 - c

        if len(mags) == prev_len:
            break  # converged

    return mags

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 65)
    print(f"  Pecaut-Mamajek Magnitude Estimator — {TARGET_NAME}")
    print("=" * 65)

    # Fetch and parse table
    rows = fetch_table(TABLE_URL)

    # Determine Teff
    teff = TEFF
    sptype_matched = None
    if teff is None and SPTYPE is not None:
        teff = sptype_to_teff(rows, SPTYPE)
        if teff is None:
            print(f"  ERROR: Spectral type '{SPTYPE}' not found in table.")
            sys.exit(1)
        print(f"  Spectral type {SPTYPE} -> Teff = {teff:.0f} K")
    elif teff is not None:
        print(f"  Input Teff = {teff:.0f} K")
    else:
        print("  ERROR: Provide either TEFF or SPTYPE in configuration.")
        sys.exit(1)

    # Find nearest spectral type
    teffs = np.array([r["Teff"] for r in rows])
    nearest_idx = np.argmin(np.abs(teffs - teff))
    nearest_row = rows[nearest_idx]
    print(f"  Nearest spectral type: {nearest_row['SpT']} "
          f"(Teff = {nearest_row['Teff']:.0f} K, "
          f"ΔTeff = {abs(nearest_row['Teff']-teff):.0f} K)")

    # Interpolate
    r1, r2 = find_nearest_rows(rows, teff)
    pm_row = interpolate_row(r1, r2, teff)

    # Check known mags
    known = {k: v for k, v in KNOWN_MAGS.items() if v is not None}
    if not known:
        print("  ERROR: No known magnitudes provided. Fill in at least one in KNOWN_MAGS.")
        sys.exit(1)
    print(f"  Known magnitudes: {known}")

    # Distance modulus
    mu, d_pc = distance_modulus(PARALLAX_MAS, DISTANCE_PC)
    if mu is not None:
        print(f"  Distance: {d_pc:.1f} pc  |  Distance modulus: {mu:.3f} mag")

    # Estimate all magnitudes
    print()
    all_mags = estimate_all_mags(known, pm_row)

    # Print results
    print("-" * 65)
    print(f"  {'Band':<8} {'Apparent mag':>13}  {'Absolute mag':>13}  {'Source'}")
    print("-" * 65)

    band_order = ["U","B","V","Rc","Ic","G","J","H","Ks","W1","W2","W3","W4"]
    band_labels = {
        "U":"Johnson U","B":"Johnson B","V":"Johnson V",
        "Rc":"Cousins Rc","Ic":"Cousins Ic","G":"Gaia G",
        "J":"2MASS J","H":"2MASS H","Ks":"2MASS Ks",
        "W1":"WISE W1","W2":"WISE W2","W3":"WISE W3","W4":"WISE W4"
    }

    for band in band_order:
        if band in all_mags:
            app_mag = all_mags[band]
            abs_mag = app_mag - mu if mu is not None else None
            source  = "INPUT" if band in known else "estimated"
            app_str = f"{app_mag:>8.3f}" if app_mag is not None else "      —"
            abs_str = f"{abs_mag:>8.3f}" if abs_mag is not None else "      —"
            label   = band_labels.get(band, band)
            print(f"  {band:<8} {app_str:>13}  {abs_str:>13}  [{source}]  {label}")

    # Print PM table color indices used
    print()
    print("-" * 65)
    print(f"  Pecaut-Mamajek interpolated parameters at Teff = {teff:.0f} K:")
    print("-" * 65)
    pm_display = [
        ("SpT",   nearest_row["SpT"],          "Nearest spectral type"),
        ("Teff",  f"{pm_row.get('Teff',np.nan):.0f} K", "Effective temperature"),
        ("Mv",    f"{pm_row.get('Mv',np.nan):.3f}",      "Abs. V magnitude"),
        ("B-V",   f"{pm_row.get('B-V',np.nan):.3f}",     "Johnson B-V"),
        ("V-Ks",  f"{pm_row.get('V-Ks',np.nan):.3f}",   "V - 2MASS Ks"),
        ("J-H",   f"{pm_row.get('J-H',np.nan):.3f}",    "2MASS J-H"),
        ("H-Ks",  f"{pm_row.get('H-Ks',np.nan):.3f}",   "2MASS H-Ks"),
        ("G-V",   f"{pm_row.get('G-V',np.nan):.3f}",    "Gaia G - V"),
        ("Ks-W1", f"{pm_row.get('Ks-W1',np.nan):.3f}",  "Ks - WISE W1"),
        ("R_Rsun",f"{pm_row.get('R_Rsun',np.nan):.3f}", "Radius (R_sun)"),
        ("Msun",  f"{pm_row.get('Msun',np.nan):.3f}",   "Mass (M_sun)"),
        ("logL",  f"{pm_row.get('logL',np.nan):.3f}",   "log(L/L_sun)"),
    ]
    for key, val, desc in pm_display:
        print(f"  {key:<10} {val:<12}  {desc}")

    print("=" * 65)
    print("  Reference: Pecaut & Mamajek (2013, ApJS, 208, 9)")
    print("  Table: https://www.pas.rochester.edu/~emamajek/")
    print("=" * 65)


if __name__ == "__main__":
    main()
