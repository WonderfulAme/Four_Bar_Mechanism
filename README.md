# Four-Bar Linkage Simulator

A standalone Python script featuring a graphical user interface (GUI) to simulate, analyze, and visualize four-bar linkages with position and velocity analysis.

---

## Features

- **Interactive GUI:** Adjust link lengths, coupler point parameters, system rotation, angular velocity, and kinematic configuration via sliders and radio buttons.
- **Mechanism Classification:** Automatically determines whether the linkage is Grashof, Non‑Grashof, or invalid, and updates status in real time.
- **Preset Mechanisms:** Load common linkage configurations (double‑crank, crank‑rocker, double‑rocker, non‑Grashof) with a single click.
- **Position Simulation:** Sweep the input crank through its valid range and record joint positions and coupler point trajectories.
- **Velocity Analysis:** Compute angular and linear velocities of all moving joints and plot velocity graphs.
- **Rotation Transformation:** Apply a global rotation to all points and velocities via a 2×2 rotation matrix.
- **Visualization:** Plot individual mechanism frames or full trajectory plots for points A, B, and P over time.

---

## Dependencies

- **Python 3.7+** (includes `tkinter` for GUI)
- **NumPy** for numerical calculations
- **Matplotlib** for plotting

Install required libraries with pip:
```bash
pip install numpy matplotlib
```

> On some Linux distributions, you may need to install `tkinter` separately:
> ```bash
> sudo apt-get install python3-tk
> ```

---

## Installation

1. **Download** `four_bar.py` to your working directory.
2. Ensure dependencies are installed (see above).

---

## Usage

### Launch the GUI

Run the script directly from the command line:
```bash
python four_bar.py
```

This opens a window with:
- **Control Panel** (left): Sliders for link lengths \(L₀\, L₁\, L₂\, L₃\), coupler length (Lap), coupler angle (α), system rotation (φ₉), angular velocity (ω₁), and configuration selection.
- **Preset Menu:** Choose a preset mechanism and load it.
- **Play/Stop Buttons:** Animate the mechanism motion.
- **Analysis View:** Switch between Mechanism animation, Position vs. time graphs, and Velocity vs. time graphs.
- **Status Display:** Shows current Grashof classification in color (green/orange/red).

### Programmatic Import

Although primarily GUI-driven, you can import and use core classes in your own scripts:

```python
import numpy as np
from four_bar import FourBarLinkage, classify_four_bar, presets

# Example: Simulate a crank-rocker
link = FourBarLinkage()
# Update parameters directly
link.update_parameters(L0=40, L1=10, L2=30, L3=50, alpha=np.radians(45), phi_g=0)
# Run simulation
link.calculate_positions()
# Access trajectories
print(link.traj_A, link.traj_B, link.traj_P)
```

---

## Script Structure

### Functions

- **`cosine_law_for_angle(adj1, adj2, opp)`**: Returns the angle opposite to `opp` using the Law of Cosines.
- **`cosine_law_for_opp_side(adj1, adj2, angle)`**: Returns the opposite side length given two adjacent sides and the included angle.
- **`classify_four_bar(lengths)`**: Standalone helper to classify a length set as "Grashof", "Non-Grashof", or "Change Point".

### `FourBarLinkage` Class

Encapsulates kinematic analysis without GUI:
- **Initialization:** Sets default link lengths, coupler settings, and computes initial positions & velocities.
- **`calculate_positions()`**: Computes joint angles and Cartesian positions for links and coupler.
- **`calculate_velocities()`**: Derives linear velocities and angular transmission as functions of time.
- **`rotate()`**: Applies a rotation matrix to all stored positions and velocities.
- **`update_parameters(**kwargs)`**: Modify any attributes (e.g., `L1`, `alpha`) and re-run calculations.
- **`plot_simulation(ax, index)`**: Draws the linkage at a single frame plus full trajectories.
- **`plot_position_graphs(fig)`**, **`plot_velocity_graphs(fig)`**: Create subplots for time-series analysis.
- **`Grashof()`**, **`is_valid_four_bar()`**, **`valid_theta_1()`**: Classification utilities.

### `FourBarLinkageSimulator` Class

Builds the Tkinter GUI around `FourBarLinkage`:
- Constructs control panels, sliders, buttons, and canvas.
- Synchronizes GUI input with kinematic model updates.
- Manages animation loops and switching between analysis views.
- Handles preset loading, error states, and status coloring.


---

## Presets

Below are the five built‑in presets, described only by their bar lengths and resulting mechanism type:

- **Grashof Double Crank**: Bar lengths L₀=12, L₁=56, L₂=67, L₃=37 → Grashof double‑crank mechanism.

- **Grashof Crank Rocker**: Bar lengths L₀=40, L₁=10, L₂=30, L₃=50 → Grashof crank‑rocker mechanism.

- **Grashof Double Rocker**: Bar lengths L₀=60, L₁=40, L₂=50, L₃=20 → Grashof double‑rocker mechanism.

- **Non‑Grashof 1**: Bar lengths L₀=55, L₁=45, L₂=35, L₃=30 → Non‑Grashof oscillating mechanism.

- **Non‑Grashof 2**: Bar lengths L₀=70, L₁=50, L₂=40, L₃=93 → Non‑Grashof oscillating mechanism.

---

## License

This project is released under the MIT License. © Amelia Hoyos Vélez

