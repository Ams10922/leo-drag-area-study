"""
Analyse VarArea_StateLog.txt from Drag_VariableArea.script.

Compares the variable-area "truth" spacecraft against the two fixed-area
models and reports the position error in the truth spacecraft's RIC frame:
  R = radial, I = in-track (along-track), C = cross-track.
In-track error is the number flight dynamics teams care about most: it is the
timing error of where the satellite will be.

Usage:
    python analyse_variable_area.py VarArea_StateLog.txt
    python analyse_variable_area.py VarArea_StateLog.txt --label "phi=0, 200 km"

Needs numpy and matplotlib (pip install numpy matplotlib).
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np

R_EARTH = 6378.1363  # km, same equatorial radius GMAT uses for JGM2

# Colours: truth = blue, fixed mean = orange, fixed max = aqua
C_VAR, C_MEAN, C_MAX = "#257de0", "#f47230", "#00b894"
INK, GRID = "#333333", "#e6e6e6"


def load(path):
    d = np.loadtxt(path)
    out = {
        "t": d[:, 0],
        "area": d[:, 1],
        "alt": d[:, 2:5],          # geodetic altitude: Var, FixMean, FixMax
        "hmean": d[:, 5:8] - R_EARTH,  # SMA - R: smoother "mean altitude"
        "var": d[:, 8:14],
        "fixmean": d[:, 14:20],
        "fixmax": d[:, 20:26],
    }
    return out


def orbit_average(t, y):
    """Running mean over ~1 orbit, to remove the J2 wobble in osculating SMA."""
    dt = np.median(np.diff(t)) * 86400.0
    a0 = np.nanmean(y[: max(1, len(y) // 100)]) + R_EARTH
    period = 2 * np.pi * np.sqrt(a0**3 / 398600.4418)
    n = max(1, int(round(period / dt)))
    out = np.full_like(y, np.nan)
    for j in range(y.shape[1]):
        col = y[:, j]
        good = ~np.isnan(col)
        c = np.convolve(np.where(good, col, 0.0), np.ones(n), "valid")
        w = np.convolve(good.astype(float), np.ones(n), "valid")
        avg = np.where(w == n, c / np.maximum(w, 1), np.nan)
        out[n // 2 : n // 2 + len(avg), j] = avg
    return out


def ric_error(truth, other):
    """Position of `other` relative to `truth`, in truth's RIC frame (km)."""
    r, v = truth[:, :3], truth[:, 3:]
    dr = other[:, :3] - r
    R = r / np.linalg.norm(r, axis=1, keepdims=True)
    h = np.cross(r, v)
    C = h / np.linalg.norm(h, axis=1, keepdims=True)
    rn = np.linalg.norm(r, axis=1)
    ro = other[:, :3]
    # In-track as an arc length along the orbit (signed angle x radius), so it
    # stays correct when the separation grows to hundreds of km or more.
    ang = np.arctan2((np.cross(r, ro) * C).sum(1), (r * ro).sum(1))
    radial = np.linalg.norm(ro, axis=1) - rn
    return np.column_stack([radial, rn * ang, (dr * C).sum(1)])


def at_day(t, y, day):
    if day > t[-1] + 1e-3:
        return np.nan
    return y[np.argmin(np.abs(t - day))]


def style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, color=INK)
    ax.set_ylabel(ylabel, color=INK)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.tick_params(colors=INK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logfile")
    ap.add_argument("--label", default="", help="text added to the figure title")
    ap.add_argument("--out", default="variable_area_comparison.png")
    ap.add_argument("--altstop", type=float, default=125.0,
                    help="AltStop used in the GMAT script (km)")
    args = ap.parse_args()

    d = load(args.logfile)
    t = d["t"]
    err_mean = ric_error(d["var"], d["fixmean"])
    err_max = ric_error(d["var"], d["fixmax"])

    # The fixed-max spacecraft is dropped from the propagation once it reaches
    # AltStop, so its state is frozen after that. Blank it out from there on.
    landed = np.nonzero(d["alt"][:, 2] <= args.altstop)[0]
    if landed.size:
        k = landed[0]
        err_max[k + 1:] = np.nan
        d["hmean"][k + 1:, 2] = np.nan
        print(f"Fixed-max spacecraft reached {args.altstop:.0f} km at day {t[k]:.3f} "
              f"(run ends at day {t[-1]:.3f})")

    # ---- Text summary --------------------------------------------------
    hbar = orbit_average(t, d["hmean"])
    last = np.nonzero(~np.isnan(hbar[:, 0]))[0][-1]
    print(f"Run length: {t[-1]:.3f} days ({len(t)} samples)")
    print(f"Area: min {d['area'].min():.3f}, max {d['area'].max():.3f}, "
          f"time-average {d['area'].mean():.4f} m^2")
    print(f"Orbit-averaged SMA - R at day {t[last]:.2f}: Var {hbar[last,0]:.2f} km, "
          f"FixMean {hbar[last,1]:.2f} km")
    dv = hbar[0 + np.nonzero(~np.isnan(hbar[:,0]))[0][0], 0] - hbar[last, 0]
    dm = hbar[0 + np.nonzero(~np.isnan(hbar[:,1]))[0][0], 1] - hbar[last, 1]
    print(f"Decay over whole run: Var {dv:.2f} km, FixMean {dm:.2f} km "
          f"({100*(dv/dm-1):+.2f}%, includes the runaway once one is lower)")
    # Clean measure of the area bias: decay over the first 1.5 days, while
    # both spacecraft are still at almost the same altitude.
    i0 = np.nonzero(~np.isnan(hbar[:, 0]))[0][0]
    i1 = np.argmin(np.abs(t - min(1.5, t[last])))
    ev = hbar[i0, 0] - hbar[i1, 0]
    em = hbar[i0, 1] - hbar[i1, 1]
    print(f"EFFECTIVE DRAG BIAS (decay to day {t[i1]:.2f}): truth decays "
          f"{100*(ev/em-1):+.2f}% vs the mean-area model "
          f"(+ = variable area behaves like a bigger area than its average)")
    print("\nIn-track error relative to the variable-area truth (km):")
    print(f"{'day':>6} {'fixed mean':>12} {'fixed max':>12}")
    for day in (0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30):
        if day <= t[-1] + 1e-3:
            print(f"{day:>6} {at_day(t, err_mean[:,1], day):>12.3f} "
                  f"{at_day(t, err_max[:,1], day):>12.3f}")
    print(f"\nMax |in-track| error: fixed mean {np.nanmax(np.abs(err_mean[:,1])):.3f} km, "
          f"fixed max {np.nanmax(np.abs(err_max[:,1])):.3f} km")
    print(f"Max |radial| error:   fixed mean {np.nanmax(np.abs(err_mean[:,0])):.3f} km, "
          f"fixed max {np.nanmax(np.abs(err_max[:,0])):.3f} km")

    # ---- Figure: four panels, one measure each -------------------------
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    title = "Variable vs fixed drag area"
    if args.label:
        title += f" ({args.label})"
    fig.suptitle(title, color=INK, fontsize=14)

    # 1. The area model itself, first ~3 orbits
    ax = axs[0, 0]
    m = t <= 0.2
    ax.plot(t[m] * 24, d["area"][m], color=C_VAR, lw=2, label="variable (truth)")
    ax.axhline(d["area"].mean(), color=C_MEAN, lw=2, ls="--", label="fixed mean")
    ax.axhline(d["area"].max(), color=C_MAX, lw=2, ls=":", label="fixed max")
    style(ax, "Elapsed time (hours)", "Drag area (m²)")
    ax.set_title("Drag area model, first few orbits", color=INK, fontsize=11)
    ax.legend(frameon=False)

    # 2. Mean altitude decay
    ax = axs[0, 1]
    ax.plot(t, hbar[:, 0], color=C_VAR, lw=2, label="variable (truth)")
    ax.plot(t, hbar[:, 1], color=C_MEAN, lw=2, ls="--", label="fixed mean")
    ax.plot(t, hbar[:, 2], color=C_MAX, lw=2, ls=":", label="fixed max")
    style(ax, "Elapsed time (days)", "Orbit-averaged SMA − R (km)")
    ax.set_title("Orbit decay", color=INK, fontsize=11)
    ax.legend(frameon=False)

    # 3. In-track error
    ax = axs[1, 0]
    ax.plot(t, err_mean[:, 1], color=C_MEAN, lw=2, label="fixed mean − truth")
    ax.plot(t, err_max[:, 1], color=C_MAX, lw=2, label="fixed max − truth")
    style(ax, "Elapsed time (days)", "In-track error (km)")
    ax.set_title("In-track (along-track) position error", color=INK, fontsize=11)
    ax.legend(frameon=False)

    # 4. In-track error of the fixed-mean model alone (its own scale)
    ax = axs[1, 1]
    ax.plot(t, err_mean[:, 1], color=C_MEAN, lw=2)
    style(ax, "Elapsed time (days)", "In-track error (km)")
    ax.set_title("Fixed mean area only (zoomed)", color=INK, fontsize=11)

    fig.savefig(args.out, dpi=150)
    print(f"\nSaved figure: {args.out}")


if __name__ == "__main__":
    main()
