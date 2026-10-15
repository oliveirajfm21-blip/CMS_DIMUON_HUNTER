"""
CMS dimuon invariant mass spectrum
==================================

Data: CMS Open Data, "Events with two muons from 2011"
      (DoubleMu 2011A, 100k events, CERN Open Data record 5201).

Steps:
    1. Configuration
    2. Load data and check integrity
    3. Compute the invariant mass (two methods)
    4. Validate against the official mass column
    5. Plot the spectrum (log-log histogram)
    6. Animated 3D event display of a random J/psi event (Plotly, HTML)

Usage:
    python main.py               # random J/psi event
    python main.py --seed 42     # reproduce a specific event
    python main.py --no-browser  # do not open the display
"""

import argparse
import webbrowser
from pathlib import Path
import urllib.request

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scienceplots  # noqa: F401  (registers the "science" plot styles)


# =========================================================
# 1. Configuration
# =========================================================

DATA_URL = "https://opendata.cern.ch/record/5201/files/Dimuon_DoubleMu.csv"
DATA_FILE = Path("data/Dimuon_DoubleMu.csv")
FIG_DIR = Path("figures")
RESULTS_DIR = Path("results")

EXPECTED_SHAPE = (100_000, 21)  # rows x columns, from the record metadata
PT_TOLERANCE = 0.01             # GeV, for the pT consistency check

# Known resonances (PDG masses, GeV) used to label the spectrum
RESONANCES = {
    r"$\rho/\omega$": 0.78,
    r"$\phi$": 1.02,
    r"$J/\psi$": 3.097,
    r"$\psi(2S)$": 3.686,
    r"$\Upsilon$": 9.46,
    r"$Z$": 91.19,
}

# Simplified CMS geometry for the event display (metres)
B_FIELD = 3.8            # T, solenoid magnetic field along +z
TRACKER_R = 1.1          # tracker outer radius
TRACKER_HALF_Z = 2.9     # tracker half-length
SOLENOID_R = 2.95        # solenoid inner radius
SOLENOID_HALF_Z = 6.25   # solenoid half-length
N_FRAMES = 60            # frames of the track animation


# ---------------------------------------------------------
# In plain words (1):
#   Before anything else, we write down every "magic number" in one
#   place: where the data lives, where the outputs go, how big the
#   file must be, the masses of the particles we expect to find and
#   the size of the detector. If something changes, we change it here
#   and nowhere else.
# ---------------------------------------------------------


# =========================================================
# 2. Load data and check integrity
# =========================================================

def load_data() -> pd.DataFrame:
    """Download the CSV if needed, read it and check it is complete."""
    if not DATA_FILE.exists():
        DATA_FILE.parent.mkdir(exist_ok=True)
        print(f"Downloading {DATA_URL} ...")
        urllib.request.urlretrieve(DATA_URL, DATA_FILE)

    df = pd.read_csv(DATA_FILE)

    # A truncated download still loads, so check against the metadata
    if df.shape != EXPECTED_SHAPE:
        raise ValueError(
            f"Incomplete file: shape {df.shape}, expected {EXPECTED_SHAPE}. "
            f"Delete {DATA_FILE} and run again."
        )
    if df.isna().sum().sum() != 0:
        raise ValueError("Missing values found: the file may be corrupted.")

    print(f"Loaded {len(df)} events, {df.shape[1]} columns.")
    return df


# ---------------------------------------------------------
# In plain words (2):
#   We open the CSV (downloading it first if it is not on disk).
#   A broken download still opens without any error, it just has
#   fewer rows. That really happened during this project: the first
#   file had only 32,592 of the 100,000 events. So the program now
#   refuses to continue unless the table has exactly the size the
#   CERN record promises and no empty cells.
# ---------------------------------------------------------


# =========================================================
# 3. Invariant mass (two methods)
# =========================================================

def mass_from_pt_eta_phi(df: pd.DataFrame) -> pd.Series:
    """M from (pT, eta, phi), assuming massless muons.

    M^2 = 2 pT1 pT2 [cosh(eta1 - eta2) - cos(phi1 - phi2)]
    """
    m2 = 2 * df["pt1"] * df["pt2"] * (
        np.cosh(df["eta1"] - df["eta2"]) - np.cos(df["phi1"] - df["phi2"])
    )
    return np.sqrt(m2)


def mass_from_four_vectors(df: pd.DataFrame) -> pd.Series:
    """Exact M from the four-vectors (E, px, py, pz), no approximation.

    M^2 = (E1 + E2)^2 - |p1 + p2|^2
    """
    m2 = (
        (df["E1"] + df["E2"]) ** 2
        - (df["px1"] + df["px2"]) ** 2
        - (df["py1"] + df["py2"]) ** 2
        - (df["pz1"] + df["pz2"]) ** 2
    )
    # Rounding can make M^2 slightly negative: mark those events as NaN
    return np.sqrt(m2.where(m2 >= 0))


# ---------------------------------------------------------
# In plain words (3):
#   Each event has two muons. If they came from the decay of one
#   heavier particle, combining their energy and momentum gives back
#   the mass of that "parent" particle. We compute it in two ways:
#     - from pT, eta and phi (the detector's natural coordinates),
#       pretending the muon has zero mass;
#     - from E, px, py, pz directly, with no approximation.
#   Both use the same idea: the energy left over after removing the
#   motion of the pair is its mass (E^2 - p^2 = M^2).
# ---------------------------------------------------------


# =========================================================
# 4. Validation against the official mass column M
# =========================================================

def validate(df: pd.DataFrame) -> str:
    """Compare both computed masses with M and check pT consistency."""
    lines = []

    # 4a. Difference to the official mass, for each method
    for col, name in [("M_calc", "pT, eta, phi (massless)"),
                      ("M_full", "four-vector (exact)")]:
        diff = (df[col] - df["M"]).abs()
        lines.append(
            f"[{name}] mean |dM| = {diff.mean():.5f} GeV, "
            f"max |dM| = {diff.max():.5f} GeV, "
            f"mean relative = {(diff / df['M']).mean():.2e}"
        )

    # 4b. Events where the pT-based mass disagrees by more than 1 GeV
    n_big = ((df["M_calc"] - df["M"]).abs() > 1).sum()
    lines.append(f"Events with |M_calc - M| > 1 GeV: {n_big}")

    # 4c. Is the pt column consistent with sqrt(px^2 + py^2)?
    for i in (1, 2):
        pt_xy = np.hypot(df[f"px{i}"], df[f"py{i}"])
        n_bad = ((pt_xy - df[f"pt{i}"]).abs() > PT_TOLERANCE).sum()
        lines.append(
            f"Muon {i}: pt inconsistent with px, py in {n_bad} events "
            f"({100 * n_bad / len(df):.1f}%)"
        )

    # 4d. Events dropped because of negative M^2
    lines.append(f"Events with negative M^2 (dropped): {df['M_full'].isna().sum()}")

    report = "\n".join(lines)
    print(report)

    RESULTS_DIR.mkdir(exist_ok=True)
    (RESULTS_DIR / "validation.txt").write_text(report + "\n", encoding="utf-8")
    return report


# ---------------------------------------------------------
# In plain words (4):
#   We check our numbers against the mass that CMS already put in the
#   file. The exact method agrees almost perfectly (only rounding).
#   The pT-based method agrees on average, but a few hundred events
#   are far off. Digging in, we found the reason: in ~3% of the muons
#   the "pt" column does not match sqrt(px^2 + py^2). So the problem
#   is in that input column, not in the formula. Everything is saved
#   to results/validation.txt so the numbers can be reproduced.
# ---------------------------------------------------------


# =========================================================
# 5. Spectrum plot (log-log histogram)
# =========================================================

def plot_spectrum(mass: pd.Series, path: Path) -> None:
    """Log-log histogram of the dimuon mass with resonance labels."""
    # Log-spaced bins: same visual width on a log axis, 0.3 to 300 GeV
    bins = np.logspace(np.log10(0.3), np.log10(300), 300)

    plt.style.use(["science", "no-latex"])  # "no-latex": no LaTeX install needed
    fig, ax = plt.subplots(figsize=(7, 4.5))

    ax.hist(mass, bins=bins, histtype="step", color="k", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Dimuon invariant mass $M_{\mu\mu}$ [GeV]")
    ax.set_ylabel("Events / bin")
    ax.set_title("CMS Open Data 2011 - DoubleMu (100k events)")

    # Vertical line and label for each known resonance
    y_top = ax.get_ylim()[1]
    for label, m in RESONANCES.items():
        ax.axvline(m, color="tab:red", linestyle=":", linewidth=0.6)
        ax.text(m, y_top * 0.5, label, rotation=90,
                ha="right", va="top", fontsize=7, color="tab:red")

    path.parent.mkdir(exist_ok=True)
    fig.savefig(path, dpi=300)
    print(f"Saved {path}")
    plt.close(fig)  # save only, do not block the script with a window


# ---------------------------------------------------------
# In plain words (5):
#   We count how many events fall in each mass interval and draw it.
#   Both axes are logarithmic because the masses go from 0.3 to 300
#   GeV and the counts from 1 to thousands. Every sharp peak is a real
#   particle that decayed into two muons (J/psi, Upsilon, Z, ...).
#   The broad bump around 10-40 GeV is not a particle: it comes from
#   the trigger, which only kept events with fast enough muons.
# ---------------------------------------------------------


# =========================================================
# 6. 3D event display (random J/psi event, animated)
# =========================================================

def kinematics_from_momentum(px, py, pz):
    """pT, eta and phi recomputed from (px, py, pz).

    Used instead of the pt/eta/phi columns, which are inconsistent
    with the momentum components in ~3% of the muons.
    """
    pt = np.hypot(px, py)        # pT = sqrt(px^2 + py^2)
    phi = np.arctan2(py, px)     # azimuthal angle, in (-pi, pi]
    eta = np.arcsinh(pz / pt)    # from pz = pT sinh(eta)
    return pt, eta, phi


def pick_jpsi_event(df: pd.DataFrame, seed: int) -> pd.Series:
    """Random J/psi candidate whose muons curve visibly inside the detector."""
    pt1, eta1, _ = kinematics_from_momentum(df["px1"], df["py1"], df["pz1"])
    pt2, eta2, _ = kinematics_from_momentum(df["px2"], df["py2"], df["pz2"])

    mask = (
        df["M_full"].between(3.05, 3.15)                # J/psi mass window
        & pt1.between(3, 8) & pt2.between(3, 8)         # moderate pT -> visible curvature
        & (eta1.abs() < 1.2) & (eta2.abs() < 1.2)       # central, inside the barrel
    )
    candidates = df[mask]
    if candidates.empty:
        raise ValueError("No J/psi event matches the display selection.")

    # Same seed -> same event, so any display can be reproduced
    event = candidates.sample(n=1, random_state=seed).iloc[0]
    print(f"Picked row {event.name} out of {len(candidates)} J/psi candidates "
          f"(seed = {seed}).")
    return event


def inside_detector(x, y, z) -> np.ndarray:
    """True while the point is inside the solenoid volume."""
    return (np.hypot(x, y) <= SOLENOID_R) & (np.abs(z) <= SOLENOID_HALF_Z)


def cut_at_exit(x, y, z):
    """Keep the trajectory points until the first exit from the solenoid."""
    inside = inside_detector(x, y, z)
    stop = len(x) if inside.all() else np.argmin(inside)
    return x[:stop], y[:stop], z[:stop]


def helix(pt, eta, phi, q, n_points=500):
    """Muon trajectory in a uniform field B along +z.

    R = pT / (0.3 B): radius in m, pT in GeV, B in T.
    s is the arc length in the transverse plane; the direction
    rotates as phi(s) = phi - q s / R (clockwise for q > 0).
    """
    R = pt / (0.3 * B_FIELD)
    s = np.linspace(0, np.pi * R, n_points)   # up to half a turn
    turn = phi - q * s / R

    x = q * R * (np.sin(phi) - np.sin(turn))
    y = q * R * (np.cos(turn) - np.cos(phi))
    z = s * np.sinh(eta)                      # dz/ds = pz / pT = sinh(eta)
    return cut_at_exit(x, y, z)


def straight_line(eta, phi, n_points=200):
    """Trajectory without magnetic field, for comparison.

    theta = 2 arctan(exp(-eta));
    direction = (sin(theta) cos(phi), sin(theta) sin(phi), cos(theta))
    """
    theta = 2 * np.arctan(np.exp(-eta))
    t = np.linspace(0, 2 * SOLENOID_HALF_Z, n_points)
    x = t * np.sin(theta) * np.cos(phi)
    y = t * np.sin(theta) * np.sin(phi)
    z = t * np.cos(theta)
    return cut_at_exit(x, y, z)


def path_length(x, y, z) -> np.ndarray:
    """Distance travelled along the trajectory up to each point (m)."""
    step = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2 + np.diff(z) ** 2)
    return np.concatenate([[0.0], np.cumsum(step)])


def build_muons(event: pd.Series) -> list:
    """Kinematics and trajectory of both muons of the event."""
    muons = []
    for i, color in ((1, "crimson"), (2, "dodgerblue")):
        pt, eta, phi = kinematics_from_momentum(
            event[f"px{i}"], event[f"py{i}"], event[f"pz{i}"]
        )
        q = int(event[f"Q{i}"])
        x, y, z = helix(pt, eta, phi, q)
        muons.append(dict(
            label="μ+" if q > 0 else "μ−", color=color,
            pt=pt, eta=eta, phi=phi, R=pt / (0.3 * B_FIELD),
            x=x, y=y, z=z, s=path_length(x, y, z),
            straight=straight_line(eta, phi),
        ))
    return muons


def event_info(event: pd.Series, muons: list, seed: int) -> str:
    """Text of the information panel (HTML line breaks)."""
    lines = [
        f"<b>Run</b> {int(event['Run'])}   <b>Event</b> {int(event['Event'])}",
        f"<b>M(μμ)</b> = {event['M_full']:.3f} GeV   (J/ψ: 3.097 GeV)",
    ]
    for m in muons:
        lines.append(
            f"<b>{m['label']}</b>: pT = {m['pt']:.2f} GeV, η = {m['eta']:.2f}, "
            f"φ = {m['phi']:.2f} rad, R = {m['R']:.2f} m"
        )
    lines.append(f"<i>seed = {seed}</i>")
    return "<br>".join(lines)


def cylinder(radius, half_z, color, opacity) -> go.Surface:
    """Transparent cylinder centred on the beam axis (z)."""
    angle = np.linspace(0, 2 * np.pi, 60)
    z = np.array([-half_z, half_z])
    angle, z = np.meshgrid(angle, z)
    return go.Surface(
        x=radius * np.cos(angle), y=radius * np.sin(angle), z=z,
        colorscale=[[0, color], [1, color]], showscale=False,
        opacity=opacity, hoverinfo="skip",
    )


def plot_event_display(event: pd.Series, path: Path, seed: int) -> None:
    """Animated 3D view: detector, collision point and both muons."""
    muons = build_muons(event)
    fig = go.Figure()

    # 6a. Static part: detector, beam axis, collision point, no-field lines
    fig.add_trace(cylinder(TRACKER_R, TRACKER_HALF_Z, "royalblue", 0.10))
    fig.add_trace(cylinder(SOLENOID_R, SOLENOID_HALF_Z, "gray", 0.06))
    fig.add_trace(go.Scatter3d(
        x=[0, 0], y=[0, 0], z=[-SOLENOID_HALF_Z, SOLENOID_HALF_Z],
        mode="lines", line=dict(color="gray", width=2, dash="dash"),
        name="Beam axis",
    ))
    fig.add_trace(go.Scatter3d(
        x=[0], y=[0], z=[0], mode="markers",
        marker=dict(size=5, color="gold"), name="Collision point",
    ))
    for m in muons:
        x, y, z = m["straight"]
        fig.add_trace(go.Scatter3d(
            x=x, y=y, z=z, mode="lines",
            line=dict(color=m["color"], width=2, dash="dot"),
            name=f"{m['label']}  without B field",
        ))

    # 6b. Animated part: each muon's track and its "head" (current position)
    first_animated = len(fig.data)
    for m in muons:
        fig.add_trace(go.Scatter3d(
            x=m["x"], y=m["y"], z=m["z"], mode="lines",
            line=dict(color=m["color"], width=6),
            name=f"{m['label']}  pT = {m['pt']:.2f} GeV",
        ))
        fig.add_trace(go.Scatter3d(
            x=m["x"][-1:], y=m["y"][-1:], z=m["z"][-1:], mode="markers",
            marker=dict(size=4, color=m["color"]), showlegend=False,
        ))
    animated = list(range(first_animated, len(fig.data)))

    # 6c. Frames: both muons move at ~ the speed of light, so in every
    #     frame they have travelled the same distance along their paths
    d_max = max(m["s"][-1] for m in muons)
    frames = []
    for k in range(1, N_FRAMES + 1):
        d = d_max * k / N_FRAMES
        data = []
        for m in muons:
            n = np.searchsorted(m["s"], d, side="right")  # points already reached
            data.append(go.Scatter3d(x=m["x"][:n], y=m["y"][:n], z=m["z"][:n]))
            data.append(go.Scatter3d(x=m["x"][n - 1:n], y=m["y"][n - 1:n],
                                     z=m["z"][n - 1:n]))
        frames.append(go.Frame(data=data, traces=animated, name=str(k)))
    fig.frames = frames

    # 6d. Layout: title, information panel and Play / Pause buttons
    play = dict(frame=dict(duration=40, redraw=True),
                transition=dict(duration=0), mode="immediate")
    pause = dict(frame=dict(duration=0, redraw=False),
                 transition=dict(duration=0), mode="immediate")
    fig.update_layout(
        title="CMS 2011 open data — J/ψ → μ+μ− event display",
        scene=dict(xaxis_title="x [m]", yaxis_title="y [m]",
                   zaxis_title="z [m] (beam)", aspectmode="data",
                   camera=dict(eye=dict(x=0.3, y=0.3, z=1.8))),  # view from above
        legend=dict(x=0, y=1),
        margin=dict(l=0, r=0, t=50, b=0),
        updatemenus=[dict(
            type="buttons", direction="left", x=0, y=0,
            xanchor="left", yanchor="bottom",
            buttons=[
                dict(label="▶ Play", method="animate", args=[None, play]),
                dict(label="❚❚ Pause", method="animate", args=[[None], pause]),
            ],
        )],
    )
    fig.add_annotation(
        text=event_info(event, muons, seed), showarrow=False,
        xref="paper", yref="paper", x=1, y=1, xanchor="right", yanchor="top",
        align="left", font=dict(size=12),
        bgcolor="rgba(255,255,255,0.85)", bordercolor="lightgray", borderwidth=1,
    )

    path.parent.mkdir(exist_ok=True)
    fig.write_html(path, include_plotlyjs="cdn")   # small file; loads Plotly from CDN
    print(f"Saved {path}")


# ---------------------------------------------------------
# In plain words (6):
#   We pick a random real J/psi event and redraw it inside a
#   simplified CMS detector. The magnetic field (3.8 T, along the
#   beam) bends each muon into a helix: slower muons (small pT) bend
#   more, and positive and negative muons bend in opposite directions.
#   The dotted lines show where each muon would go with no field.
#   Pressing Play makes both muons leave the collision point at the
#   same speed (close to the speed of light), so one may leave the
#   detector before the other. A panel shows the event's numbers.
#   The random choice uses a "seed": the same seed always gives the
#   same event, so any picture can be reproduced exactly.
# ---------------------------------------------------------


# =========================================================
# Main
# =========================================================

def main(seed=None, open_browser=True) -> None:
    df = load_data()

    df["M_calc"] = mass_from_pt_eta_phi(df)
    df["M_full"] = mass_from_four_vectors(df)

    validate(df)

    # The four-vector mass is the reference (the pt column has inconsistencies)
    plot_spectrum(df["M_full"].dropna(), FIG_DIR / "dimuon_spectrum.png")

    # No seed given: draw one, and print it so the event can be reproduced
    if seed is None:
        seed = int(np.random.default_rng().integers(1_000_000))
    display_path = FIG_DIR / "event_display_jpsi.html"
    plot_event_display(pick_jpsi_event(df, seed), display_path, seed)

    if open_browser:
        webbrowser.open(display_path.resolve().as_uri())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CMS dimuon spectrum and J/psi event display")
    parser.add_argument("--seed", type=int, default=None,
                        help="seed for the displayed event (default: random)")
    parser.add_argument("--no-browser", action="store_true",
                        help="do not open the event display in the browser")
    args = parser.parse_args()
    main(seed=args.seed, open_browser=not args.no_browser)


# ---------------------------------------------------------
# In plain words (main):
#   main() runs the steps in order: load -> compute the mass ->
#   validate -> plot the spectrum -> draw one random event and open
#   it in the browser. The line `if __name__ == "__main__"` makes this
#   happen only when the file is run directly (python main.py), not
#   when another script or a notebook imports one of its functions.
#   argparse reads the optional --seed and --no-browser flags.
# ---------------------------------------------------------
