# nirps-mag-estimator

Estimate the photometric magnitudes of a main-sequence star in any band,
starting from its effective temperature (or spectral type) and **at least one**
known magnitude.

The script interpolates the colour indices of the
[Pecaut & Mamajek (2013, ApJS, 208, 9)](https://ui.adsabs.harvard.edu/abs/2013ApJS..208....9P/abstract)
dwarf colour/T<sub>eff</sub> sequence and propagates them from the bands you
provide to every other band in the table.

Written for target selection and exposure-time estimates for
[NIRPS](https://www.eso.org/sci/facilities/lasilla/instruments/nirps.html)
observing programmes, where the H-band magnitude of a candidate is often
missing, but it is generic and works for any main-sequence star.

## Bands supported

| System | Bands |
|---|---|
| Johnson–Cousins | U, B, V, R<sub>c</sub>, I<sub>c</sub> |
| Gaia | G |
| 2MASS | J, H, K<sub>s</sub> |
| WISE | W1, W2, W3, W4 |

If a parallax or distance is given, absolute magnitudes are reported alongside
the apparent ones.

## Installation

```bash
git clone https://github.com/<your-username>/nirps-mag-estimator.git
cd nirps-mag-estimator
pip install -r requirements.txt
```

Requires Python 3.8+ and `numpy`. The Pecaut–Mamajek table is downloaded from
[Eric Mamajek's page](https://www.pas.rochester.edu/~emamajek/) on first run;
if you have no internet access, download
`EEM_dwarf_UBVIJHK_colors_Teff.txt` manually and place it in the working
directory — the script will pick it up automatically.

## Usage

Edit the `CONFIGURATION` block at the top of `nirps_mag_estimator.py`, then run:

```bash
python nirps_mag_estimator.py
```

```python
TARGET_NAME  = "TOI-4349"

TEFF         = 3850      # K — takes priority over SPTYPE
SPTYPE       = None      # e.g. "K2V", "M3V" — used when TEFF is None

KNOWN_MAGS = {           # fill in what you have, leave the rest as None
    "J": 10.4,
}

PARALLAX_MAS = 13.89     # or set DISTANCE_PC instead (optional)
DISTANCE_PC  = None
```

### Example output

```
=================================================================
  Pecaut-Mamajek Magnitude Estimator — TOI-4349
=================================================================
  Input Teff = 3850 K
  Nearest spectral type: M0V (Teff = 3870 K, ΔTeff = 20 K)
  Known magnitudes: {'J': 10.4}
  Distance: 72.0 pc  |  Distance modulus: 4.287 mag

-----------------------------------------------------------------
  Band     Apparent mag   Absolute mag  Source
-----------------------------------------------------------------
  V              13.442         9.155   [estimated]  Johnson V
  G              12.610         8.323   [estimated]  Gaia G
  J              10.400         6.113   [INPUT]      2MASS J
  H               9.760         5.473   [estimated]  2MASS H
  Ks              9.575         5.288   [estimated]  2MASS Ks
  ...
```

(Numbers shown are illustrative — run the script for the real values.)

## Caveats

- Valid for **main-sequence (dwarf) stars only**. Giants, subgiants and
  pre-main-sequence stars follow different colour sequences.
- No extinction or reddening correction is applied. For reddened targets the
  derived colours — and therefore the estimated magnitudes —will be biased.
- Colours are linearly interpolated between the two nearest table rows in
  T<sub>eff</sub>; accuracy degrades at the extreme ends of the table.
- Uncertainties are not propagated. Treat the output as an estimate good to
  roughly a few tenths of a magnitude, not a measurement.

## Contributing

Issues and pull requests are welcome. Useful things to add:

- a command-line interface, so the config block need not be edited
- uncertainty propagation from the input magnitude and T<sub>eff</sub>
- optional extinction correction
- a vectorised mode for whole target lists

## Citation

If this is useful in published work, please cite the underlying table:

> Pecaut, M. J. & Mamajek, E. E. 2013, ApJS, 208, 9

## License

MIT — see [LICENSE](LICENSE).
