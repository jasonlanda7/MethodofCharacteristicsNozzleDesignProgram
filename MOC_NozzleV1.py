# moc_bell_nozzle.py
# =========================================================================
#  METHOD OF CHARACTERISTICS (MOC) FOR AXISYMMETRIC NOZZLE DESIGN (Python)
# =========================================================================
#
#  Author: [Jason Da Silva]  |  Python port by: ChatGPT
#  Credits: Characteristic-line formulation & base logic adapted from
#           VDEngineering on YouTube (2025).
#
#  Date: 2025-11-12
#  Version: 1.0 (Python)
#
#  PURPOSE:
#  --------
#  Designs a shock-free, minimum-length, axisymmetric supersonic nozzle
#  using the Method of Characteristics (MOC) to construct the internal
#  expansion field and resulting nozzle wall contour from the throat (M=1)
#  to a target exit Mach number implied by P0/Pe (Pe set below).
#
#  Inputs:
#     P0 : stagnation pressure (Pa)
#     g  : ratio of specific heats (gamma)
#     TR : throat radius (same units as you want in output plots; e.g. mm)
#
#  Returns (dict):
#     {
#       "xw": wall x (upper) up to exit,
#       "yw": wall y (upper) up to exit,
#       "xw_mirror": mirrored x (lower),
#       "yw_mirror": mirrored y (lower),
#       "full_x": concatenated full contour x (upper + mirrored),
#       "full_y": concatenated full contour y,
#       "Me": computed exit Mach number,
#       "T_max_deg": maximum turning angle in degrees
#     }
#
#  Notes:
#   - Replicates the MATLAB algorithm closely (angles in rad internally).
#   - Writes two sheets to PARAMS.xlsx ("PTS" and "FullContour"),
#     matching the MATLAB writematrix calls.
# =========================================================================

from __future__ import annotations

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar


def moc_bell_nozzle(P0: float, g: float, TR: float, do_plot: bool = True):
    # Ambient / exit pressure (design condition)
    Pe = 101325.0  # Pa

    # Exit Mach from isentropic P0/Pe relation
    Me = math.sqrt((2.0 / (g - 1.0)) * ((P0 / Pe) ** ((g - 1.0) / g) - 1.0))

    # Degrees <-> radians
    DtoR = math.pi / 180.0
    RtoD = 180.0 / math.pi

    # Prandtl–Meyer function ν(M) for perfect gas (radians)
    A = math.sqrt((g + 1.0) / (g - 1.0))
    B = (g - 1.0) / (g + 1.0)

    def V_PM(M):
        # guard for M>1
        if M <= 1.0:
            return 0.0
        return A * math.atan(math.sqrt(B * (M * M - 1.0))) - math.atan(math.sqrt(M * M - 1.0))

    # Max wall turning angle (radians) for minimum-length nozzle
    T_max = 0.5 * V_PM(Me)  # radians
    T_max_deg = T_max * RtoD

    # MATLAB splits the first ray angle using the fractional part DT of (90 - T_max_deg)
    DT = (90.0 - T_max_deg) - math.floor(90.0 - T_max_deg)
    # first theta entry in radians
    T = [DT * DtoR]

    # number of interior rays n = 2*T_max_deg (as in MATLAB)
    n = int(round(T_max_deg * 2.0))

    # storage
    P = [None]  # 1-index placeholder mimic; we’ll drop later to match MATLAB indexing
    M_list = [None]
    RR = [None]
    LR = [None]
    SL = [None]

    # Build the fan (rays) and their slopes
    # T(m) in radians, m = 2..n+1 (MATLAB style)
    for m in range(2, n + 2):
        Tm = (DT + (m - 1)) * DtoR   # radians
        T.append(Tm)

        # Solve T(m) - V_PM(M) = 0 for M in (1, 1.01*Me) via bracketing
        # This mirrors MATLAB's fzero([1, 1.01*Me])
        def f(M):
            return Tm - V_PM(M)

        lo = 1.0 + 1e-8
        hi = max(1.01 * Me, 1.0001)  # ensure hi>lo

        # If bracket fails, back off conservatively by widening upper bound
        bracket_found = False
        for scale in [1.01, 1.05, 1.1, 1.2, 1.5, 2.0]:
            try:
                res = root_scalar(f, bracket=[lo, hi * scale], method="brentq")
                M_sol = float(res.root)
                bracket_found = True
                break
            except ValueError:
                pass
        if not bracket_found:
            # fallback: Newton with a reasonable guess
            M_guess = max(1.001, Me * 0.9)
            res = root_scalar(f, x0=M_guess, fprime=None, method="newton")
            M_sol = float(res.root)

        M_list.append(M_sol)

        # X-axis point where C+ from centerline hits wall (MATLAB: P(m) = 0 + TR*tan(T(m)))
        Pm = TR * math.tan(Tm)
        P.append(Pm)

        # Right-running (C+) slope from centerline point to wall point: RR = -TR/P
        RR.append(-TR / Pm)

        # Left-running (C−) slope from grid point: LR = tan(T + μ), μ = asin(1/M)
        mu = math.asin(1.0 / M_sol)
        LR.append(math.tan(Tm + mu))

        # SL = -RR (used later)
        SL.append(-RR[-1])

    # Drop the initial placeholder to match MATLAB’s post-processing (P(1) = [])
    P = P[1:]
    RR = RR[1:]
    LR = LR[1:]
    SL = SL[1:]
    T = T  # first entry is T[0] = DT*DtoR, same as MATLAB T(1)

    # Construct intersections of interior characteristics with the initial fan
    # F = RR(end) in MATLAB (used as a constant slope)
    F = RR[-1]

    # Intersections (x, y) for interior nodes along last C+ slope F
    x = []
    y = []
    for c in range(len(P) - 1):
        # x(c) = (TR + SL(c)*P(c)) / (SL(c) - F)
        xc = (TR + SL[c] * P[c]) / (SL[c] - F)
        yc = F * xc + TR
        x.append(xc)
        y.append(yc)

    # === Wall Section Generation ===
    TM = T_max  # radians
    # First wall intersection (with the first C− through P[0])
    # xw(1) = (TR + SL(1)*P(1)) / (SL(1) - tan(TM))
    xw = []
    yw = []
    xw1 = (TR + SL[0] * P[0]) / (SL[0] - math.tan(TM))
    yw1 = math.tan(TM) * xw1 + TR
    xw.append(xw1)
    yw.append(yw1)

    # Initialize wall slope ladder: s(k) decreases linearly from tan(TM)
    DTW = math.tan(TM) / (len(P) - 1)
    s = [math.tan(TM)]
    b = [TR]

    # March the wall points
    for k in range(1, len(P) - 1):  # k = 2..length(P)-1 in MATLAB; zero-based here
        sk = math.tan(TM) - k * DTW
        s.append(sk)

        # b(k) = yw(k-1) - s(k)*xw(k-1)
        bk = yw[-1] - sk * xw[-1]
        b.append(bk)

        # xw(k) = (b(k) + SL(k)*P(k)) / (SL(k) - s(k))
        xwk = (bk + SL[k] * P[k]) / (SL[k] - sk)
        ywk = sk * xwk + bk

        xw.append(xwk)
        yw.append(ywk)

    # Second pass in MATLAB re-plots lines from interior nodes to wall points.
    # Geometry already captured above; replicate end state only.

    # Determine nozzle exit from final characteristic
    k_end = len(P) - 1
    x_exit = xw[k_end - 1]  # last valid index in zero-based (k_end-1)
    y_exit = yw[k_end - 1]

    # Trim/append (matching MATLAB’s behavior)
    xw = xw[:k_end]
    yw = yw[:k_end]
    xw.append(x_exit)
    yw.append(y_exit)

    # Mirror for axisymmetric visualization
    xw_mirror = xw.copy()
    yw_mirror = (-np.array(yw)).tolist()

    # Full contour (upper then mirrored back)
    full_x = xw + list(reversed(xw_mirror))
    full_y = yw + list(reversed(yw_mirror))

    # === Excel output ===
    # "PTS" sheet: columns A(xw), B(yw)
    # "FullContour" sheet: columns A(full_x), B(full_y)
    with pd.ExcelWriter("PARAMS.xlsx", engine="openpyxl") as xlw:
        # === Wall points (upper half only)
        pts_df = pd.DataFrame({
            "Axial Distance x [mm]": xw,
            "Radius y [mm]": yw
        })
        pts_df.to_excel(xlw, sheet_name="PTS", index=False)

        # === Full mirrored contour
        full_df = pd.DataFrame({
            "Axial Distance x [mm]": full_x,
            "Radius y [mm]": full_y
        })
        full_df.to_excel(xlw, sheet_name="FullContour", index=False)

    # === Plotting (optional) ===
    if do_plot:
        plt.figure(figsize=(8, 6))
        # Initial fan lines (centerline to P[j])
        for Pj in P:
            plt.plot([Pj, 0.0], [0.0, TR], linewidth=1.0)

        # Intersections to wall (interior grid to wall)
        for j in range(len(x)):
            plt.plot([x[j], xw[min(j + 1, len(xw) - 1)]],
                     [y[j], yw[min(j + 1, len(yw) - 1)]],
                     linewidth=1.0)

        # Upper wall and mirrored wall
        plt.plot(xw, yw, linewidth=2.0)
        plt.plot(xw_mirror, yw_mirror, linewidth=2.0)

        # Exit and centerline
        plt.plot([x_exit, x_exit], [y_exit, -y_exit], linestyle="--", linewidth=1.5)
        plt.axhline(0.0, linestyle="--", linewidth=1.0)

        plt.gca().set_aspect("equal", adjustable="box")
        plt.grid(True)
        plt.xlabel("Axial Distance (same units as TR)")
        plt.ylabel("Radius (same units as TR)")
        plt.title("Full Axisymmetric Nozzle Contour (Method of Characteristics)")
        plt.tight_layout()
        plt.show()

    # --- Derived geometric quantities ---
    exit_radius = yw[-1]
    nozzle_length = xw[-1]

    return {
        "xw": np.array(xw),
        "yw": np.array(yw),
        "xw_mirror": np.array(xw_mirror),
        "yw_mirror": np.array(yw_mirror),
        "full_x": np.array(full_x),
        "full_y": np.array(full_y),
        "Me": Me,
        "T_max_deg": T_max_deg,
        "exit_radius": exit_radius,
        "nozzle_length": nozzle_length,
    }


if __name__ == "__main__":
    print("\n=== AXISYMMETRIC NOZZLE DESIGN (MOC) ===")

    # Prompt user for inputs
    P0 = float(input("Enter chamber (stagnation) pressure P0 [Pa]: "))
    gamma = float(input("Enter specific heat ratio (gamma): "))
    TR = float(input("Enter throat radius [mm]: "))

    # Optional plotting flag
    plot_choice = input("Display plot? (y/n): ").strip().lower()
    do_plot = plot_choice in ['y', 'yes']

    # Call the solver
    results = moc_bell_nozzle(P0, gamma, TR, do_plot)

    # Display summary
    # Display summary
    print("\n--- Results ---")
    print(f"Exit Mach number (Me): {results['Me']:.4f}")
    print(f"Max turning angle (deg): {results['T_max_deg']:.3f}")
    print(f"Exit radius (Re): {results['exit_radius']:.3f} mm")
    print(f"Nozzle length (L): {results['nozzle_length']:.3f} mm")
    print("Data written to: PARAMS.xlsx\n")
