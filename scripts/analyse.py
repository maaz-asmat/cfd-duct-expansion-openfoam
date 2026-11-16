"""Post-process all cases: friction validation, loss coefficients, mesh independence, plots."""
import sys, pathlib, csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from case_builder import DuctCase
from foam_io import axial_profile, line_fit, read_internal_field, latest_time
from run_cases import CASES

OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
ref = DuctCase()
H1, H2, U1, nu = ref.H1, ref.H2, ref.U1, ref.nu
sigma = H1 / H2                                  # area ratio
U2 = U1 * sigma
q1 = 0.5 * U1 ** 2                               # kinematic dynamic pressure, inlet duct
Dh1, Dh2 = 2 * H1, 2 * H2
RHO = 998.0


def dean_f(U, H):
    """Darcy friction factor of a smooth wide duct from Dean's correlation Cf = 0.073 Re_H^-0.25."""
    return 4 * 0.073 * (U * H / nu) ** -0.25


K_BC = (1 - sigma) ** 2                          # Borda-Carnot, uniform profiles
K_BC_corr = 1.04 - 2 * 1.015 * sigma + (2 * 1.015 - 1.04) * sigma ** 2   # with turbulent alpha, beta


def analyse(name):
    c = CASES[name]
    d = ROOT / "cases" / name
    info = json.loads((d / "done.json").read_text())
    x, p0, ps = axial_profile(d)
    Ld = c.diffuser_length
    up = line_fit(x, p0, -6 * H1, -1 * H1)
    dn = line_fit(x, p0, Ld + 25 * H2, Ld + 57 * H2)
    K = (np.polyval(up, 0.0) - np.polyval(dn, Ld)) / q1
    ups = line_fit(x, ps, -6 * H1, -1 * H1)
    dns = line_fit(x, ps, Ld + 25 * H2, Ld + 57 * H2)
    Cp = (np.polyval(dns, Ld) - np.polyval(ups, 0.0)) / q1
    t = latest_time(d)
    C = read_internal_field(t / "C"); U = read_internal_field(t / "U")
    near = (C[:, 1] > 0.97 * H2 / 2) & (C[:, 0] > Ld) & (C[:, 0] < Ld + 2.0)
    back = near & (U[:, 0] < 0)
    L_rec = float(C[back, 0].max() - Ld) if back.any() else 0.0
    # convergence: change of mean inlet pressure over the last 500 iterations
    hist = sorted((ROOT / "cases" / name / "postProcessing").rglob("surfaceFieldValue*.dat"))
    vals = np.vstack([np.loadtxt(h, comments="#") for h in hist])
    vals = vals[vals[:, 0] <= int(t.name)]
    last = vals[-1, 1]; prev = vals[vals[:, 0] <= vals[-1, 0] - 500][-1, 1]
    return dict(case=name, angle=c.diffuser_angle_deg, mesh_scale=c.mesh_scale, cells=info["cells"],
                iterations=int(t.name), K=K, Cp=Cp, f_up=-up[0] * Dh1 / q1, f_dn=-dn[0] * Dh2 / (0.5 * U2 ** 2),
                L_diffuser_mm=1000 * Ld, L_recirc_mm=1000 * L_rec,
                conv_change_pct=100 * abs(last - prev) / q1, profile=(x, p0, ps, up, dn))


res = {n: analyse(n) for n in CASES if (ROOT / "cases" / n / "done.json").exists()}
for r in res.values():
    print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items() if k != "profile"})

# ---------- tables
s = res["sudden"]
with open(OUT / "01_validation.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["quantity", "CFD", "reference", "reference_source", "difference_%"])
    w.writerow(["friction factor, inlet duct (Re_Dh = 100 000)", round(s["f_up"], 5), round(dean_f(U1, H1), 5), "Dean (1978)", round(100 * (s["f_up"] / dean_f(U1, H1) - 1), 1)])
    w.writerow(["friction factor, outlet duct (Re_Dh = 100 000)", round(s["f_dn"], 5), round(dean_f(U2, H2), 5), "Dean (1978)", round(100 * (s["f_dn"] / dean_f(U2, H2) - 1), 1)])
    w.writerow(["sudden-expansion loss K", round(s["K"], 4), round(K_BC, 4), "Borda-Carnot, uniform profiles", round(100 * (s["K"] / K_BC - 1), 1)])
    w.writerow(["sudden-expansion loss K", round(s["K"], 4), round(K_BC_corr, 4), "Borda-Carnot with turbulent profile corrections", round(100 * (s["K"] / K_BC_corr - 1), 1)])

mesh = [res[n] for n in ["sudden_coarse", "sudden_medium", "sudden"] if n in res]
with open(OUT / "02_mesh_independence.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["case", "mesh_scale", "cells", "K", "friction_factor_inlet", "recirculation_length_mm"])
    for r in mesh:
        w.writerow([r["case"], r["mesh_scale"], r["cells"], round(r["K"], 4), round(r["f_up"], 5), round(r["L_recirc_mm"])])

order = ["sudden", "diffuser_60deg", "diffuser_30deg", "diffuser_20deg", "diffuser_10deg"]
var = [res[n] for n in order if n in res]
Ks = s["K"]
with open(OUT / "03_design_variants.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["variant", "total_angle_deg", "diffuser_length_mm", "K", "K_vs_sudden_%", "pressure_recovery_Cp",
                "ideal_Cp", "loss_Pa_at_1m_per_s_water", "recirculation_length_mm", "inlet_p_change_last500_iter_%_of_dyn_pressure"])
    for r in var:
        w.writerow([r["case"], r["angle"], round(r["L_diffuser_mm"]), round(r["K"], 3), round(100 * (r["K"] / Ks - 1), 1),
                    round(r["Cp"], 3), round(1 - sigma ** 2, 3), round(r["K"] * q1 * RHO), round(r["L_recirc_mm"]), round(r["conv_change_pct"], 3)])

# ---------- plot 1: method
x, p0, ps, up, dn = s["profile"]
fig, ax = plt.subplots(figsize=(9, 4.3))
m = (x > -0.7) & (x < 5.9)
ax.plot(x[m], RHO * p0[m], color="#1f5f8b", label="total pressure (mass-flow weighted)")
ax.plot(x[m], RHO * ps[m], color="#d97706", label="static pressure (area weighted)")
xx = np.linspace(-0.7, 0.4, 10); ax.plot(xx, RHO * np.polyval(up, xx), "k--", lw=1, label="fully developed friction lines")
xx = np.linspace(-0.4, 5.9, 10); ax.plot(xx, RHO * np.polyval(dn, xx), "k--", lw=1)
ax.annotate("", (0.03, RHO * np.polyval(dn, 0)), (0.03, RHO * np.polyval(up, 0)), arrowprops=dict(arrowstyle="<->", color="crimson"))
ax.text(0.12, RHO * (np.polyval(up, 0) + np.polyval(dn, 0)) / 2, f"expansion loss\nK = {Ks:.3f}", color="crimson", va="center")
ax.axvline(0, color="0.6", lw=0.8)
ax.set(xlabel="x [m]  (expansion at x = 0)", ylabel="pressure [Pa]", title="Sudden expansion, duct height 50 → 100 mm, water, Re = 100 000")
ax.grid(alpha=.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "01_pressure_along_duct.png", dpi=150)

# ---------- plot 2: loss vs angle
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ang = [r["angle"] for r in var]; K = [r["K"] for r in var]; L = [r["L_diffuser_mm"] for r in var]
ax[0].plot(ang, K, "o-", color="#0e7c7b", ms=8)
ax[0].axhline(K_BC, color="0.4", ls="--", lw=1); ax[0].text(100, K_BC + 0.006, "Borda–Carnot (sudden)", fontsize=8, color="0.3")
for a, k in zip(ang, K):
    ax[0].annotate(f"{k:.3f}", (a, k), (7, 2), textcoords="offset points", fontsize=8)
ax[0].set(xlabel="Total diffuser angle [deg]  (180 = sudden)", ylabel="Loss coefficient K (ref. inlet velocity)", title="Pressure loss vs diffuser angle")
ax[0].set_xticks([10, 20, 30, 60, 90, 120, 180]); ax[0].grid(alpha=.3)
ax[1].plot(L, K, "s-", color="#d97706", ms=7)
for a, l, k in zip(ang, L, K):
    ax[1].annotate(f"{a:.0f}°", (l, k), (6, 4), textcoords="offset points", fontsize=8)
ax[1].set(xlabel="Diffuser length [mm]", ylabel="Loss coefficient K", title="Trade-off: loss vs installation length"); ax[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(OUT / "03_loss_vs_angle.png", dpi=150)

# ---------- plot 3: velocity fields (mirrored about the centreline)
cases_plot = [(n, t) for n, t in [("sudden", "Sudden expansion"), ("diffuser_30deg", "30° diffuser"), ("diffuser_10deg", "10° diffuser")] if n in res]
fig, axs = plt.subplots(len(cases_plot), 1, figsize=(10, 2.0 * len(cases_plot) + 0.4), sharex=True)
axs = np.atleast_1d(axs)
for a, (n, title) in zip(axs, cases_plot):
    c = CASES[n]
    t = latest_time(ROOT / "cases" / n)
    C = read_internal_field(t / "C"); U = read_internal_field(t / "U")
    m = (C[:, 0] > -0.1) & (C[:, 0] < 1.0)
    xs, ys, ux = np.r_[C[m, 0], C[m, 0]], np.r_[C[m, 1], -C[m, 1]], np.r_[U[m, 0], U[m, 0]]
    tri = mtri.Triangulation(xs, ys)
    xc, yc = xs[tri.triangles].mean(axis=1), np.abs(ys[tri.triangles]).mean(axis=1)
    Ld = c.diffuser_length
    wall = np.where(xc < 0, H1 / 2, np.where(xc < Ld, H1 / 2 + (xc / max(Ld, 1e-9)) * (H2 - H1) / 2, H2 / 2))
    tri.set_mask(yc > wall * 0.995)
    cf = a.tricontourf(tri, ux, levels=np.linspace(-0.25, 1.25, 31), cmap="RdBu_r", extend="both")
    a.tricontour(tri, ux, levels=[0], colors="k", linewidths=0.8)
    a.set_ylim(-H2 / 2, H2 / 2); a.set_ylabel("y [m]")
    a.set_title(f"{title}: axial velocity (black line = 0, reverse flow inside)", fontsize=9)
axs[-1].set_xlabel("x [m]  (vertical scale stretched for visibility)")
fig.colorbar(cf, ax=axs, label="Ux [m/s]", shrink=0.9)
fig.savefig(OUT / "02_velocity_fields.png", dpi=150, bbox_inches="tight")
