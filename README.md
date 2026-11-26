# Variable vs fixed drag area in LEO orbit propagation

How much along-track accuracy do you lose by modelling a satellite as a fixed-area "cannonball" when its real cross-section changes through the orbit?

**Status:** work in progress, v0.1 (26 Nov 2026). Two array phases run at 200 km so far; see [Next steps](#next-steps).

## Key result (200 km circular orbit, i = 30°)

A fixed drag area equal to the orbit average predicted re-entry time to within about 1%, but its along-track position error reached about 100 km by day 5 and about 420 km by re-entry.

| | Reached 125 km | Along-track error at day 1 / 3 / 5 |
|---|---|---|
| Variable area (truth) | day 6.88 | n/a |
| Fixed mean area, 1.64 m² | about 1.5 h later* | 3 / 30 / 104 km behind |
| Fixed max area, 2.0 m² | day 5.70 | 64 / 721 / 2,911 km ahead |

\*From a separate fixed-area run with the same ballistic coefficient (382 kg/m²).

The variable-area spacecraft also decayed about 1% faster than its own average area predicts. The explanation is that with phi = 0 the array is face-on at the equator, where Earth's equatorial bulge puts the spacecraft about 5 km lower, in denser air. A phase test confirmed it:

| Array phase | Array face-on at | Effective drag bias | Fixed-mean along-track error, day 3 |
|---|---|---|---|
| phi = 0 | equator crossings | +1.04% | 30.0 km behind |
| phi = 90 deg | +/-30 deg latitude | -1.01% | 29.6 km ahead |

Shifting where the area peaks reverses the bias and flips the error, so the timing of the area relative to the density, not just its average, sets the effective drag.

![Variable vs fixed drag area, 200 km, phi = 0](results/200km_phi0.png)

The full write-up is in the technical note: [docs/technical_note_v0.1.pdf](docs/technical_note_v0.1.pdf).

## Method in brief

- **Tool:** GMAT R2026a. Runge–Kutta 8(9), JGM-2 gravity at degree and order 8, NRLMSISE-00 with constant F10.7 = 150 and Kp = 3, spherical drag model, Cd = 2.2.
- **Three spacecraft**, identical except for drag area, propagated together from the same state:
  - **Truth:** area recomputed every 60 s as `A = A_body + A_panel * |cos(u + phi)|`, a body plus a rotating array (u = argument of latitude).
  - **Fixed mean:** the orbit-average area, `A_body + A_panel * 2/pi`.
  - **Fixed max:** the array always face-on, `A_body + A_panel`.
- **Error metric:** in-track error is the signed arc length along the truth's orbit. Decay is the semi-major axis minus Earth radius, averaged over one orbit to remove the J2 oscillation.

## Repository layout

```
gmat/Drag_VariableArea.script      GMAT mission: three spacecraft, variable-area loop (phi = 0, 10 days)
gmat/Drag_VariableArea_phi90.script  Same, phi = 90 deg, 3 days
analysis/analyse_variable_area.py  Reads the GMAT log, prints the summary, plots the figure
results/                           Figures from each run
docs/                              Technical note (PDF)
```

## How to run

1. Open `gmat/Drag_VariableArea.script` in GMAT (tested with R2026a on Windows).
2. Edit the variables at the top of the mission sequence if needed: `Abody`, `Apanel`, `phi`, `Tend`, `AltStop`. To change altitude, set `SMA` on all three spacecraft to `6378.1363 + h`.
3. Run it. GMAT writes `VarArea_StateLog.txt` to its `output` folder; a 7-day run at 200 km takes about 7 minutes.
4. Analyse the log (needs Python 3 with numpy and matplotlib):

   ```
   python analysis/analyse_variable_area.py VarArea_StateLog.txt --label "200 km, phi = 0" --out results/200km_phi0.png
   ```

   This prints the in-track errors and the "effective drag bias" (how much faster the truth decays than the fixed-mean model over the first 1.5 days), and saves a four-panel figure.

## Next steps

- [ ] Null test: `Apanel = 0`, all three orbits must coincide
- [x] Phase test: phi = 90 deg (bias -1.0%, as predicted)
- [ ] phi = 45 and 135 deg (a difference between them would show the day-night density effect)
- [ ] 500 km case over 30 days
- [ ] Area error vs density error: F10.7 = 70 and 250
- [ ] Attitude-based area table (box + panels, nadir pointing) using GMAT's SPAD drag model
- [ ] Cross-check one case in Orekit

## Limitations

This is a toy area model (one body, one array, a |cos| profile) with a fixed drag coefficient and constant space weather, and so far only one altitude, inclination and phase. See section 4 of the technical note.

## Background

The question came out of very-low-Earth-orbit work for the Airbus Space Competition 2026/27, where an earlier ballistic-coefficient sweep had used a fixed spherical drag area.
