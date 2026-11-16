"""
Builds an OpenFOAM (v2412, simpleFoam, k-omega SST) case for turbulent water flow
through a 2D duct expansion: channel height H1 = 50 mm -> H2 = 100 mm (area ratio 0.5).

Only the upper half of the symmetric duct is modelled (symmetry plane on the centreline).
Expansion types:
  - sudden expansion        (diffuser_angle_deg = 180)
  - planar diffuser         (total included angle 2*theta, e.g. 10, 20, 30, 60 deg)

The inlet profile is recycled from 8*H1 downstream, so the flow reaching the
expansion is fully developed. Units: SI; pressure is kinematic (p/rho).
"""
from dataclasses import dataclass
import math
import pathlib
import textwrap


@dataclass
class DuctCase:
    H1: float = 0.050               # m, inlet channel height (full)
    H2: float = 0.100               # m, outlet channel height (full)
    U1: float = 1.0                 # m/s, bulk inlet velocity
    nu: float = 1.0e-6              # m^2/s, water at ~20 C
    diffuser_angle_deg: float = 180.0   # total included angle; 180 = sudden expansion
    L_up_H1: float = 15.0           # upstream length in H1
    map_H1: float = 8.0             # inlet profile recycled from this many H1 downstream
    L_down_H2: float = 60.0         # downstream length in H2
    n_y1: int = 40                  # cells across the inlet half-channel
    n_y2: int = 40                  # extra cells across the step (sudden) / added height (diffuser)
    wall_cell: float = 2.0e-5       # m, first cell at walls (y+ ~ 1, wall-resolved SST)
    dx: float = 0.004               # m, axial cell size in the diffuser
    end_iterations: int = 3000
    mesh_scale: float = 1.0         # multiplies all cell counts (mesh-independence study)

    @property
    def Dh1(self):                  # hydraulic diameter of a wide duct = 2 x height
        return 2 * self.H1

    @property
    def Re(self):
        return self.U1 * self.Dh1 / self.nu

    @property
    def diffuser_length(self):
        if self.diffuser_angle_deg >= 180:
            return 0.0
        return (self.H2 - self.H1) / 2 / math.tan(math.radians(self.diffuser_angle_deg / 2))


HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def _w(path, cls, obj, body):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER.format(cls=cls, obj=obj) + textwrap.dedent(body))


def grading_ratio(total, n, end_cell):
    """Expansion ratio (last/first) of n geometric cells over `total` whose last cell is `end_cell`."""
    lo, hi = 1e-6, 1.0
    for _ in range(200):
        q = 0.5 * (lo + hi)
        first = end_cell / q ** (n - 1)
        length = first * (1 - q ** n) / (1 - q)
        lo, hi = (q, hi) if length > total else (lo, q)
    return q ** (n - 1)


def block_mesh(c: DuctCase) -> str:
    m = c.mesh_scale
    h1, h2 = c.H1 / 2, c.H2 / 2                       # half heights (model y-extent)
    x0 = -c.L_up_H1 * c.H1
    Ld = c.diffuser_length
    x2 = Ld + c.L_down_H2 * c.H2
    z = 0.001
    sudden = c.diffuser_angle_deg >= 180
    wc = c.wall_cell / m
    ny1, ny2 = int(round(c.n_y1 * m)), int(round(c.n_y2 * m))
    n_up, n_down = int(round(120 * m)), int(round(400 * m))

    names, pts = {}, []

    def P(name, x, y):
        names[name] = len(pts)
        pts.append((x, y))

    P("c0", x0, 0); P("w0", x0, h1); P("c1", 0, 0); P("w1", 0, h1)
    if sudden:
        P("s1", 0, h2); P("c2", x2, 0); P("w2", x2, h1); P("s2", x2, h2)
    else:
        P("c1d", Ld, 0); P("w1d", Ld, h2); P("c2", x2, 0); P("w2", x2, h2)

    nP = len(pts)

    def v(n, back=False):
        return names[n] + (nP if back else 0)

    def hexb(a, b, cc, d, cells, grading):
        ids = [v(a), v(cc), v(d), v(b), v(a, 1), v(cc, 1), v(d, 1), v(b, 1)]
        return f"    hex ({' '.join(map(str, ids))}) ({cells[0]} {cells[1]} 1) simpleGrading ({grading})\n"

    g1 = grading_ratio(h1, ny1, wc)
    blocks = ""
    if sudden:
        g_half = grading_ratio((h2 - h1) / 2, ny2 // 2, wc)
        g_step = f"((0.5 0.5 {1 / g_half:.6g}) (0.5 0.5 {g_half:.6g}))"
        blocks += hexb("c0", "w0", "c1", "w1", (n_up, ny1), f"0.3 {g1:.6g} 1")
        blocks += hexb("c1", "w1", "c2", "w2", (n_down, ny1), f"20 {g1:.6g} 1")
        blocks += hexb("w1", "s1", "w2", "s2", (n_down, ny2), f"20 {g_step} 1")
    else:
        g = grading_ratio(h1, ny1 + ny2, wc)          # same ratio along the duct so block faces match
        n_d = max(20, int(round(Ld / c.dx * m)))
        blocks += hexb("c0", "w0", "c1", "w1", (n_up, ny1 + ny2), f"0.3 {g:.6g} 1")
        blocks += hexb("c1", "w1", "c1d", "w1d", (n_d, ny1 + ny2), f"1 {g:.6g} 1")
        blocks += hexb("c1d", "w1d", "c2", "w2", (n_down, ny1 + ny2), f"20 {g:.6g} 1")

    def quad(a, b):            # face between two front points a,b and their back copies
        return f"({v(a)} {v(b)} {v(b, 1)} {v(a, 1)})"

    if sudden:
        inlet, outlet = [quad("c0", "w0")], [quad("c2", "w2"), quad("w2", "s2")]
        walls = [quad("w0", "w1"), quad("w1", "s1"), quad("s1", "s2")]
        sym = [quad("c0", "c1"), quad("c1", "c2")]
        empties = [("c0", "c1", "w1", "w0"), ("c1", "c2", "w2", "w1"), ("w1", "w2", "s2", "s1")]
    else:
        inlet, outlet = [quad("c0", "w0")], [quad("c2", "w2")]
        walls = [quad("w0", "w1"), quad("w1", "w1d"), quad("w1d", "w2")]
        sym = [quad("c0", "c1"), quad("c1", "c1d"), quad("c1d", "c2")]
        empties = [("c0", "c1", "w1", "w0"), ("c1", "c1d", "w1d", "w1"), ("c1d", "c2", "w2", "w1d")]
    fb = [f"({' '.join(str(v(n)) for n in e)})" for e in empties] + \
         [f"({' '.join(str(v(n, 1)) for n in e)})" for e in empties]

    def patch(name, typ, faces, extra=""):
        return (f"    {name}\n    {{\n        type {typ};\n{extra}        faces\n        (\n"
                + "".join(f"            {f}\n" for f in faces) + "        );\n    }\n")

    mapped = ("        sampleMode nearestCell;\n        sampleRegion region0;\n        samplePatch none;\n"
              f"        offset ({c.map_H1 * c.H1:.6g} 0 0);\n")
    vtx = "".join(f"    ({x:.9g} {y:.9g} 0)\n" for x, y in pts) + "".join(f"    ({x:.9g} {y:.9g} {z})\n" for x, y in pts)
    return (f"\nscale 1;\n\nvertices\n(\n{vtx});\n\nblocks\n(\n{blocks});\n\nedges ();\n\nboundary\n(\n"
            + patch("inlet", "mappedPatch", inlet, mapped) + patch("outlet", "patch", outlet)
            + patch("walls", "wall", walls) + patch("centreline", "symmetryPlane", sym)
            + patch("frontAndBack", "empty", fb) + ");\n")


def write_case(c: DuctCase, case_dir):
    case_dir = pathlib.Path(case_dir)
    k_in = 1.5 * (0.05 * c.U1) ** 2
    omega_in = math.sqrt(k_in) / (0.09 ** 0.25 * 0.07 * c.Dh1)

    _w(case_dir / "system/blockMeshDict", "dictionary", "blockMeshDict", block_mesh(c))
    _w(case_dir / "system/controlDict", "dictionary", "controlDict", f"""
    application     simpleFoam;
    startFrom       latestTime;
    startTime       0;
    stopAt          endTime;
    endTime         {c.end_iterations};
    deltaT          1;
    writeControl    timeStep;
    writeInterval   {c.end_iterations};
    purgeWrite      1;
    writeFormat     ascii;
    writePrecision  8;
    writeCompression off;
    timeFormat      general;
    timePrecision   6;
    runTimeModifiable true;

    functions
    {{
        inletPressure
        {{
            type            surfaceFieldValue;
            libs            (fieldFunctionObjects);
            writeControl    timeStep;
            writeInterval   50;
            writeFields     false;
            regionType      patch;
            name            inlet;
            operation       areaAverage;
            fields          (p);
        }}
    }}
    """)
    _w(case_dir / "system/fvSchemes", "dictionary", "fvSchemes", """
    ddtSchemes      { default steadyState; }
    gradSchemes     { default Gauss linear; }
    divSchemes
    {
        default         none;
        div(phi,U)      bounded Gauss linearUpwind grad(U);
        div(phi,k)      bounded Gauss limitedLinear 1;
        div(phi,omega)  bounded Gauss limitedLinear 1;
        div((nuEff*dev2(T(grad(U))))) Gauss linear;
    }
    laplacianSchemes { default Gauss linear corrected; }
    interpolationSchemes { default linear; }
    snGradSchemes   { default corrected; }
    wallDist        { method meshWave; }
    """)
    _w(case_dir / "system/fvSolution", "dictionary", "fvSolution", """
    solvers
    {
        p { solver GAMG; tolerance 1e-08; relTol 0.01; smoother GaussSeidel; }
        "(U|k|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-08; relTol 0.1; }
    }
    SIMPLE
    {
        nNonOrthogonalCorrectors 0;
        consistent      yes;
        residualControl { p 1e-6; U 1e-7; "(k|omega)" 1e-7; }
    }
    relaxationFactors { equations { U 0.9; ".*" 0.9; } }
    """)
    _w(case_dir / "constant/transportProperties", "dictionary", "transportProperties", f"""
    transportModel  Newtonian;
    nu              {c.nu};
    """)
    _w(case_dir / "constant/turbulenceProperties", "dictionary", "turbulenceProperties", """
    simulationType  RAS;
    RAS { RASModel kOmegaSST; turbulence on; printCoeffs off; }
    """)
    common = """
        centreline   { type symmetryPlane; }
        frontAndBack { type empty; }"""
    _w(case_dir / "0/U", "volVectorField", "U", f"""
    dimensions      [0 1 -1 0 0 0 0];
    internalField   uniform ({c.U1} 0 0);
    boundaryField
    {{
        inlet  {{ type mapped; field U; setAverage true; average ({c.U1} 0 0); interpolationScheme cell; value uniform ({c.U1} 0 0); }}
        outlet {{ type inletOutlet; inletValue uniform (0 0 0); value uniform ({c.U1} 0 0); }}
        walls  {{ type noSlip; }}{common}
    }}
    """)
    _w(case_dir / "0/p", "volScalarField", "p", f"""
    dimensions      [0 2 -2 0 0 0 0];
    internalField   uniform 0;
    boundaryField
    {{
        inlet  {{ type zeroGradient; }}
        outlet {{ type fixedValue; value uniform 0; }}
        walls  {{ type zeroGradient; }}{common}
    }}
    """)
    for name, val, dims, wall in [
            ("k", k_in, "[0 2 -2 0 0 0 0]", "type fixedValue; value uniform 1e-12;"),
            ("omega", omega_in, "[0 0 -1 0 0 0 0]", f"type omegaWallFunction; value uniform {omega_in:.6g};")]:
        _w(case_dir / f"0/{name}", "volScalarField", name, f"""
    dimensions      {dims};
    internalField   uniform {val:.6g};
    boundaryField
    {{
        inlet  {{ type mapped; field {name}; setAverage false; average 0; interpolationScheme cell; value uniform {val:.6g}; }}
        outlet {{ type inletOutlet; inletValue uniform {val:.6g}; value uniform {val:.6g}; }}
        walls  {{ {wall} }}{common}
    }}
    """)
    _w(case_dir / "0/nut", "volScalarField", "nut", f"""
    dimensions      [0 2 -1 0 0 0 0];
    internalField   uniform 0;
    boundaryField
    {{
        inlet  {{ type calculated; value uniform 0; }}
        outlet {{ type calculated; value uniform 0; }}
        walls  {{ type nutLowReWallFunction; value uniform 0; }}{common}
    }}
    """)
    return case_dir
