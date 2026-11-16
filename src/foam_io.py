"""Minimal readers for OpenFOAM ASCII field files and pressure-loss post-processing."""
import pathlib
import re
import numpy as np


def read_internal_field(path):
    """Return the internalField of an ASCII OpenFOAM field as a numpy array."""
    txt = pathlib.Path(path).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", txt)
    if not m:
        u = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
        raise ValueError(f"uniform field in {path}: {u.group(1) if u else '?'}")
    kind, n = m.group(1), int(m.group(2))
    body = txt[m.end():]
    if kind == "scalar":
        vals = np.array(body[:body.index(")")].split(), float)
    else:
        vals = np.array(re.findall(r"\(([^()]*)\)", body[:body.index(")\n;")])[:n] if ")\n;" in body else
                        re.findall(r"\(([^()]*)\)", body)[:n])
        vals = np.array([v.split() for v in vals], float)
    assert len(vals) == n, (path, len(vals), n)
    return vals


def latest_time(case):
    times = [p for p in pathlib.Path(case).iterdir() if p.is_dir() and re.fullmatch(r"[0-9.eE+-]+", p.name) and p.name != "0"]
    return max(times, key=lambda p: float(p.name))


def axial_profile(case, bin_width=0.002):
    """Mass-flow-weighted total pressure p0 = p + |U|^2/2 (kinematic) and area-weighted static pressure along x.
    In the 2D duct the cell volume is proportional to its cross-section area, so V is the area weight."""
    t = latest_time(case)
    C = read_internal_field(t / "C")
    V = read_internal_field(t / "V")
    U = read_internal_field(t / "U")
    p = read_internal_field(t / "p")
    x = C[:, 0]
    p0 = p + 0.5 * np.sum(U ** 2, axis=1)
    w = np.clip(U[:, 0], 0, None) * V                     # mass-flux weight (forward flow)
    wa = V                                                # area weight for static pressure
    edges = np.arange(x.min(), x.max() + bin_width, bin_width)
    k = np.digitize(x, edges)
    xs, P0, PS = [], [], []
    for i in np.unique(k):
        s = k == i
        if w[s].sum() <= 0:
            continue
        xs.append(np.average(x[s], weights=wa[s]))
        P0.append(np.average(p0[s], weights=w[s]))
        PS.append(np.average(p[s], weights=wa[s]))
    return np.array(xs), np.array(P0), np.array(PS)


def line_fit(x, y, lo, hi):
    s = (x >= lo) & (x <= hi)
    return np.polyfit(x[s], y[s], 1)
