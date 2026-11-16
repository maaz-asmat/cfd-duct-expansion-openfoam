# Repeating the study in ANSYS Fluent (Student)

The model is 2D and small, so it runs on an ordinary laptop with the free ANSYS Student licence.

## 1. Geometry (SpaceClaim or DesignModeler, 2D)
Model the upper half of the duct. All dimensions are in mm, and the expansion is at x = 0.
- Inlet duct: x from −750 to 0, y from 0 to 25.
- Sudden expansion: outlet duct from x = 0 to 6000, y from 0 to 50.
- Diffuser variants: a straight wall from (0, 25) to (L, 50), with L = 25 / tan(angle/2). That gives L = 286, 142, 93 and 43 mm for 10°, 20°, 30° and 60° total angle.

## 2. Mesh
- Mapped (structured) quad mesh. First cell at the walls **0.02 mm**, so that y⁺ ≈ 1. Use about 40 cells across the inlet half-height and grow them towards the centre.
- Use fine axial cells (2–5 mm) near the expansion and coarser ones further away.
- Named selections: `inlet`, `outlet`, `walls`, `centreline`.

## 3. Setup
- Solver: pressure-based, steady, 2D planar.
- Model: **k-omega SST**. Material: water, ρ = 998 kg/m³, μ = 0.001 kg/(m·s).
- Inlet: velocity inlet, 1 m/s, turbulence intensity 5 %, hydraulic diameter 0.1 m. For a fully developed profile, either make the inlet duct about 40 heights long, or use a periodic section first and write its profile.
- Outlet: pressure outlet, 0 Pa. Centreline: symmetry. Walls: no slip.
- Methods: Coupled, second-order upwind for momentum, k and ω.

## 4. Post-processing
1. Create vertical lines every 0.1 m along the duct. Report the **mass-weighted average total pressure** on each line and export to CSV.
2. Fit straight lines to the fully developed regions: x = −0.30 to −0.05 m upstream, and x = 2.5 to 5.7 m downstream.
3. K = (gap between the two lines at the expansion) / (0.5 × 998 × 1²).
4. Compare with `results/03_design_variants.csv`. Expect K ≈ 0.27 for the sudden expansion.
