import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
np.set_printoptions(precision=16)  # Set decimal precision to 16 places


def cosine_law_for_angle(adj1, adj2, opp):
    return np.arccos((adj1**2 + adj2**2 - opp**2) / (2 * adj1 * adj2))

def cosine_law_for_opp_side(adj1, adj2, angle):
    return np.sqrt(adj1**2 + adj2**2 - 2 * adj1 * adj2 * np.cos(angle))

class FourBarLinkage:
    def __init__(self):
        self.L0 = 50  # Ground link (fixed)
        self.L1 = 30  # Input link
        self.L2 = 10  # Coupler link
        self.L3 = 40  # Output link
        self.Lap = 5  # Length A to P
        self.alpha = np.radians(30)
        self.config = 1
        self.phi_g = np.radians(140) # Rotation angle of the system
        self.omega_1 = 30.16  # Angular velocity of input link (rad/s)
        self.calculate_positions()
        self.calculate_velocities()
        self.rotate()
        self.valid = True

    def update_parameters(self, **kwargs):
        """Update parameters and recalculate positions"""
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
        self.calculate_positions()
        self.calculate_velocities()
        self.rotate()

    def Grashof(self):
        """Check if the mechanism is Grashof, Non-Grashof"""
        lengths = [self.L0, self.L1, self.L2, self.L3]
        S = min(lengths)
        L = max(lengths)
        P, Q = sorted(lengths)[1:3]
        if not self.is_valid_four_bar():
            self.valid = False
            return "Invalid"
        elif S + L > P + Q:
            return "Non-Grashof"
        else:
            return "Grashof" 
        
    def is_valid_four_bar(self):
        lengths = [self.L0, self.L1, self.L2, self.L3]
        for i in range(4):
            if lengths[i] >= sum(lengths) - lengths[i]:
                return False
        return True


    def valid_theta_1(self, config=1):
        """Calculate valid ranges for input angle θ₁"""
        if self.Grashof() == "Invalid":
            raise ValueError("Invalid four-bar linkage configuration")
        if self.Grashof() == "Non-Grashof":
            # Non-Grashof mechanism
            if max(self.L0, self.L1, self.L2, self.L3) in [self.L0, self.L1]:
                phi = cosine_law_for_angle(self.L1, self.L0, self.L2 + self.L3)
                if 0 < phi:
                    return (-phi, phi)
                else:
                    return (-phi + 2*np.pi, phi + 2*np.pi)
            else:
                phi = cosine_law_for_angle(self.L1, self.L0, abs(self.L2 - self.L3))
                if 0 > phi:
                    return (phi - 2*np.pi, -phi)
                else:
                    return (phi, -phi + 2*np.pi)
        elif self.Grashof() == "Grashof":
            if config == 1:
                if min(self.L0, self.L1, self.L2, self.L3) in [self.L0, self.L1]:
                    return (0, 2 * np.pi)  # Full rotation possible
                else:
                    phi_1 = cosine_law_for_angle(self.L1, self.L0, abs(self.L2 - self.L3))
                    phi_2 = cosine_law_for_angle(self.L1, self.L0, abs(self.L2 + self.L3))
                    if 0 > phi_1:
                        return (-phi_2, -phi_1)
                    elif - 0 > phi_2:
                        return (-phi_2 + 2*np.pi,  - phi_1 + 2*np.pi)
                    else:
                        return (phi_1, phi_2)
            else:
                if min(self.L0, self.L1, self.L2, self.L3) in [self.L0, self.L1]:
                    return (0, 2 * np.pi)  # Full rotation possible
                else:
                    phi_1 = cosine_law_for_angle(self.L1, self.L0, abs(self.L2 - self.L3))
                    phi_2 = cosine_law_for_angle(self.L1, self.L0, abs(self.L2 + self.L3))
                    if 0 > phi_1:
                        return (phi_1, phi_2)
                    elif 0 > phi_2:
                        return (phi_1 + 2*np.pi, phi_2 + 2*np.pi)
                    else:
                        return (- phi_2 + 2*np.pi, - phi_1 + 2*np.pi)
                    
    def calculate_positions(self):
        """Calculate all positions based on current parameters."""
        # Initialize empty lists
        self.positions = []
        self.traj_A = []
        self.traj_B = []
        self.traj_P = []
        self.theta_1_list = []  # Store theta_1 values for analysis
        self.theta_2_list = []  # Store theta_2 values for analysis

        try:
            self.theta_1_valid = self.valid_theta_1(config=self.config)
            if not isinstance(self.theta_1_valid, tuple) or len(self.theta_1_valid) != 2:
                raise ValueError("Invalid theta_1 range")

            # Define theta_1 array
            self.theta_1_arr = np.linspace(self.theta_1_valid[0], self.theta_1_valid[1], 100)
            for theta_1 in self.theta_1_arr:
                    # Use cosine law to calculate linkage positions
                    L_ac = cosine_law_for_opp_side(self.L1, self.L0, theta_1)
                    beta = cosine_law_for_angle(self.L0, L_ac, self.L1)
                    psi = cosine_law_for_angle(self.L2, L_ac, self.L3)
                    lmd = cosine_law_for_angle(L_ac, self.L3, self.L2)
                    
                    theta_2 = psi - beta
                    theta_3 = np.pi - lmd - beta

                    if theta_1 > np.pi:
                        theta_2 = psi + beta
                        theta_3 = np.pi - lmd + beta


                    # Store angles for analysis
                    self.theta_1_list.append(theta_1)
                    self.theta_2_list.append(theta_2)

                    # Compute positions of the linkage points
                    Ax = self.L1 * np.cos(theta_1)
                    Ay = self.L1 * np.sin(theta_1)

                    Bx = self.L3 * np.cos(theta_3)
                    By = self.L3 * np.sin(theta_3)

                    Px = Ax + self.Lap * np.cos(theta_2 + self.alpha)
                    Py = Ay + self.Lap * np.sin(theta_2 + self.alpha)

                    Cx = self.L0
                    Cy = 0

                    # Store positions and trajectories
                    self.positions.append(((0, 0), (Ax, Ay), (Bx, By), (Cx, Cy), (Px, Py)))
                    self.traj_A.append((Ax, Ay))
                    self.traj_B.append((Bx, By))
                    self.traj_P.append((Px, Py))
        
        except Exception as e:
            print(f"Error in position calculation: {e}")
            
    def calculate_velocities(self):
        """Calculate velocities for points A, B, and P"""
        if not hasattr(self, 'theta_1_list') or not self.theta_1_list:
            return
            
        self.vel_A = []
        self.vel_B = []
        self.vel_P = []
        self.omega_2_list = []  # Angular velocity of link 2
        
        # Calculate the time vector (assuming constant angular velocity omega_1)
        self.time = []
        t = 0
        prev_theta_1 = self.theta_1_list[0]
        
        for i, theta_1 in enumerate(self.theta_1_list):
            # Handle angle wrapping for time calculation
            if i > 0:
                delta_theta = theta_1 - prev_theta_1
                # Fix for angle wrapping
                if delta_theta > np.pi:
                    delta_theta -= 2*np.pi
                elif delta_theta < -np.pi:
                    delta_theta += 2*np.pi
                t += abs(delta_theta) / self.omega_1
            
            self.time.append(t)
            prev_theta_1 = theta_1
            
            # Calculate velocity of point A
            vAx = -self.L1 * self.omega_1 * np.sin(theta_1)
            vAy = self.L1 * self.omega_1 * np.cos(theta_1)
            self.vel_A.append((vAx, vAy))
            
            # Calculate omega_2 (angular velocity of link 2)
            # First calculate the transmission angle and link positions
            theta_2 = self.theta_2_list[i]
            
            # Calculate omega_2 using velocity constraint equations
            # vBx = vAx - omega_2 * L2 * sin(theta_2)
            # vBy = vAy + omega_2 * L2 * cos(theta_2)
            # Solving for omega_2 using dot product with normal vector to link 2
            
            # Normal vector to link 2 (perpendicular to the link)
            nx = -np.sin(theta_2)
            ny = np.cos(theta_2)
            
            # Dot product of vA with normal vector divided by L2 gives omega_2
            omega_2 = (vAx * nx + vAy * ny) / self.L2
            self.omega_2_list.append(omega_2)
            
            # Calculate velocity of point B
            vBx = vAx - self.L2 * omega_2 * np.sin(theta_2)
            vBy = vAy + self.L2 * omega_2 * np.cos(theta_2)
            self.vel_B.append((vBx, vBy))
            
            # Calculate velocity of coupler point P
            # Use the same omega_2 for the coupler
            vPx = vAx - self.Lap * omega_2 * np.sin(theta_2 + self.alpha)
            vPy = vAy + self.Lap * omega_2 * np.cos(theta_2 + self.alpha)
            self.vel_P.append((vPx, vPy))

    def rotate(self):
        # Rotate the system
        # Create the rotation matrix
        R = np.array([
            [np.cos(self.phi_g), -np.sin(self.phi_g)],
            [np.sin(self.phi_g),  np.cos(self.phi_g)]
        ])
        
        # Rotate positions
        for i in range(len(self.positions)):
            position = self.positions[i]
            (O, A, B, C, P) = position
            # Define point vectors
            O = np.array([0, 0])
            A = np.array([A[0], A[1]])
            B = np.array([B[0], B[1]])
            P = np.array([P[0], P[1]])
            C = np.array([C[0], C[1]])

            # Rotate using matrix multiplication
            O = R @ O
            A = R @ A
            B = R @ B
            P = R @ P
            C = R @ C

            # Store the rotated positions
            self.positions[i] = (O, A, B, C, P)
            self.traj_A[i] = (A[0], A[1])
            self.traj_B[i] = (B[0], B[1])
            self.traj_P[i] = (P[0], P[1])
        
        # Rotate velocities
        if hasattr(self, 'vel_A') and self.vel_A:
            for i in range(len(self.vel_A)):
                vA = np.array([self.vel_A[i][0], self.vel_A[i][1]])
                vB = np.array([self.vel_B[i][0], self.vel_B[i][1]])
                vP = np.array([self.vel_P[i][0], self.vel_P[i][1]])
                
                vA = R @ vA
                vB = R @ vB
                vP = R @ vP
                
                self.vel_A[i] = (vA[0], vA[1])
                self.vel_B[i] = (vB[0], vB[1])
                self.vel_P[i] = (vP[0], vP[1])

    def plot_simulation(self, ax, index):
        """Plot the mechanism at a specific position index"""
        (O, A, B, C, P) = self.positions[index]
        ax.plot([O[0], A[0]], [O[1], A[1]], 'ro-', label="Link 1")
        ax.plot([A[0], B[0]], [A[1], B[1]], 'go-', label="Link 2")
        ax.plot([B[0], C[0]], [B[1], C[1]], 'bo-', label="Link 3")
        ax.plot([O[0], C[0]], [O[1], C[1]], 'mo-', label="Ground Link L0")
        ax.plot([A[0], P[0]], [A[1], P[1]], 'co-', label="PA")
        ax.plot(*zip(*self.traj_A), 'r--', label="Trajectory A")
        ax.plot(*zip(*self.traj_B), 'g--', label="Trajectory B")
        ax.plot(*zip(*self.traj_P), 'k--', label="Trajectory P")

        return ax

    def plot_position_graphs(self, fig):
        """Plot position graphs for A, B, and P"""
        if not hasattr(self, 'time') or not self.time:
            return
            
        # Extract x and y coordinates from trajectories
        Ax = [p[0] for p in self.traj_A]
        Ay = [p[1] for p in self.traj_A]
        Bx = [p[0] for p in self.traj_B]
        By = [p[1] for p in self.traj_B]
        Px = [p[0] for p in self.traj_P]
        Py = [p[1] for p in self.traj_P]
        
        # Create 2x2 subplot
        axs = fig.subplots(2, 2)
        
        # X position vs time
        axs[0, 0].plot(self.time, Ax, 'r-', label='A_x')
        axs[0, 0].plot(self.time, Bx, 'g-', label='B_x')
        axs[0, 0].plot(self.time, Px, 'b-', label='P_x')
        axs[0, 0].set_title('X Position vs Time')
        axs[0, 0].set_xlabel('Time (s)')
        axs[0, 0].set_ylabel('Position X (units)')
        axs[0, 0].grid(True)
        axs[0, 0].legend()
        
        # Y position vs time
        axs[0, 1].plot(self.time, Ay, 'r-', label='A_y')
        axs[0, 1].plot(self.time, By, 'g-', label='B_y')
        axs[0, 1].plot(self.time, Py, 'b-', label='P_y')
        axs[0, 1].set_title('Y Position vs Time')
        axs[0, 1].set_xlabel('Time (s)')
        axs[0, 1].set_ylabel('Position Y (units)')
        axs[0, 1].grid(True)
        axs[0, 1].legend()
        
        # Angular positions vs time
        axs[1, 0].plot(self.time, self.theta_1_list, 'r-', label='θ₁')
        axs[1, 0].plot(self.time, self.theta_2_list, 'g-', label='θ₂')
        axs[1, 0].set_title(f'Angular Positions vs Time [{np.degrees(self.theta_1_valid[0]):.2f}° to {np.degrees(self.theta_1_valid[1]):.2f}°]')
        axs[1, 0].set_xlabel('Time (s)')
        axs[1, 0].set_ylabel('Angle (rad)')
        axs[1, 0].grid(True)
        axs[1, 0].legend()
        
        # Coupler point trajectory
        axs[1, 1].plot(Px, Py, 'b-', label='P')
        axs[1, 1].set_title('Coupler Point P Trajectory')
        axs[1, 1].set_xlabel('X Position')
        axs[1, 1].set_ylabel('Y Position')
        axs[1, 1].grid(True)
        axs[1, 1].set_aspect('equal')
        axs[1, 1].legend()
        
        fig.tight_layout()
        
    def plot_velocity_graphs(self, fig):
        """Plot velocity graphs for A, B, and P"""
        if not hasattr(self, 'time') or not self.time:
            return
            
        # Extract velocity components
        vAx = [v[0] for v in self.vel_A]
        vAy = [v[1] for v in self.vel_A]
        vBx = [v[0] for v in self.vel_B]
        vBy = [v[1] for v in self.vel_B]
        vPx = [v[0] for v in self.vel_P]
        vPy = [v[1] for v in self.vel_P]
        
        # Calculate velocity magnitudes
        vA_mag = [np.sqrt(vx**2 + vy**2) for vx, vy in zip(vAx, vAy)]
        vB_mag = [np.sqrt(vx**2 + vy**2) for vx, vy in zip(vBx, vBy)]
        vP_mag = [np.sqrt(vx**2 + vy**2) for vx, vy in zip(vPx, vPy)]
        
        # Create 2x2 subplot
        axs = fig.subplots(2, 2)
        
        # X velocity vs time
        axs[0, 0].plot(self.time, vAx, 'r-', label='vA_x')
        axs[0, 0].plot(self.time, vBx, 'g-', label='vB_x')
        axs[0, 0].plot(self.time, vPx, 'b-', label='vP_x')
        axs[0, 0].set_title('X Velocity vs Time')
        axs[0, 0].set_xlabel('Time (s)')
        axs[0, 0].set_ylabel('Velocity X (units/s)')
        axs[0, 0].grid(True)
        axs[0, 0].legend()
        
        # Y velocity vs time
        axs[0, 1].plot(self.time, vAy, 'r-', label='vA_y')
        axs[0, 1].plot(self.time, vBy, 'g-', label='vB_y')
        axs[0, 1].plot(self.time, vPy, 'b-', label='vP_y')
        axs[0, 1].set_title('Y Velocity vs Time')
        axs[0, 1].set_xlabel('Time (s)')
        axs[0, 1].set_ylabel('Velocity Y (units/s)')
        axs[0, 1].grid(True)
        axs[0, 1].legend()
        
        # Velocity magnitude vs time
        axs[1, 0].plot(self.time, vA_mag, 'r-', label='|vA|')
        axs[1, 0].plot(self.time, vB_mag, 'g-', label='|vB|')
        axs[1, 0].plot(self.time, vP_mag, 'b-', label='|vP|')
        axs[1, 0].set_title('Velocity Magnitude vs Time')
        axs[1, 0].set_xlabel('Time (s)')
        axs[1, 0].set_ylabel('Velocity (units/s)')
        axs[1, 0].grid(True)
        axs[1, 0].legend()
        
        # Angular velocity vs time
        axs[1, 1].plot(self.time, [self.omega_1] * len(self.time), 'r-', label='ω₁')
        axs[1, 1].plot(self.time, self.omega_2_list, 'g-', label='ω₂')
        axs[1, 1].set_title('Angular Velocity vs Time')
        axs[1, 1].set_xlabel('Time (s)')
        axs[1, 1].set_ylabel('Angular Velocity (rad/s)')
        axs[1, 1].grid(True)
        axs[1, 1].legend()
        
        fig.tight_layout()
            
class FourBarLinkageSimulator:
    def __init__(self, root):
        self.root = root
        self.root.title("Four-Bar Linkage Simulator with Position and Velocity Analysis")
        self.root.geometry("1200x800")
        
        # Initialize default parameters
        self.L0 = 10  # Ground link (fixed)
        self.L1 = 20  # Input link
        self.L2 = 30  # Coupler link
        self.L3 = 40  # Output link
        self.Lap = 5  # Length A to P
        self.alpha = np.radians(30)  # Angle for point P
        self.phi_g = np.radians(140)  # System rotation angle
        self.config = 1  # Configuration (1 or -1)
        self.omega_1 = 1.0  # Angular velocity of input link (rad/s)
        self.animation_running = False
        self.animation_speed = 50
        self.current_frame = 0
        self.analysis_mode = "Mechanism"  # Default view
        
        # Create frames
        self.control_frame = tk.Frame(root)
        self.control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        
        self.plot_frame = tk.Frame(root)
        self.plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Create matplotlib figures
        self.fig = Figure(figsize=(8, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Create controls
        self.create_controls()
        
        # Initialize linkage and update display
        self.update_linkage()
        
    def create_controls(self):
        # Title
        title_label = tk.Label(self.control_frame, text="Four-Bar Linkage", font=("Arial", 12, "bold"))
        title_label.grid(row=0, column=0, pady=10, sticky="w")

        # Info Button
        info_btn = tk.Button(self.control_frame, text="ℹ️ License", command= show_license_info)
        info_btn.grid(row=0, column=1, pady=10, sticky="e")

        # Link length controls
        tk.Label(self.control_frame, text="Ground Link (L₀):").grid(row=1, column=0, sticky=tk.W, pady=3)
        self.L0_slider = tk.Scale(self.control_frame, from_=10, to=100, orient=tk.HORIZONTAL, length=200,
                                 command=lambda x: self.update_parameter('L0', float(x)))
        self.L0_slider.set(self.L0)
        self.L0_slider.grid(row=1, column=1, pady=3)
        
        tk.Label(self.control_frame, text="Input Link (L₁):").grid(row=2, column=0, sticky=tk.W, pady=3)
        self.L1_slider = tk.Scale(self.control_frame, from_=10, to=100, orient=tk.HORIZONTAL, length=200,
                                 command=lambda x: self.update_parameter('L1', float(x)))
        self.L1_slider.set(self.L1)
        self.L1_slider.grid(row=2, column=1, pady=3)
        
        tk.Label(self.control_frame, text="Coupler Link (L₂):").grid(row=3, column=0, sticky=tk.W, pady=3)
        self.L2_slider = tk.Scale(self.control_frame, from_=10, to=100, orient=tk.HORIZONTAL, length=200,
                                 command=lambda x: self.update_parameter('L2', float(x)))
        self.L2_slider.set(self.L2)
        self.L2_slider.grid(row=3, column=1, pady=5)
        
        tk.Label(self.control_frame, text="Output Link (L₃):").grid(row=4, column=0, sticky=tk.W, pady=3)
        self.L3_slider = tk.Scale(self.control_frame, from_=10, to=100, orient=tk.HORIZONTAL, length=200,
                                 command=lambda x: self.update_parameter('L3', float(x)))
        self.L3_slider.set(self.L3)
        self.L3_slider.grid(row=4, column=1, pady=3)
        
        # Point P controls
        tk.Label(self.control_frame, text="Point P Length:").grid(row=5, column=0, sticky=tk.W, pady=3)
        self.Lap_slider = tk.Scale(self.control_frame, from_=0, to=50, orient=tk.HORIZONTAL, length=200,
                                  command=lambda x: self.update_parameter('Lap', float(x)))
        self.Lap_slider.set(self.Lap)
        self.Lap_slider.grid(row=5, column=1, pady=3)
        
        tk.Label(self.control_frame, text="Point P Angle (°):").grid(row=6, column=0, sticky=tk.W, pady=3)
        self.alpha_slider = tk.Scale(self.control_frame, from_=0, to=360, orient=tk.HORIZONTAL, length=200,
                                   command=lambda x: self.update_parameter('alpha', np.radians(float(x))))
        self.alpha_slider.set(np.degrees(self.alpha))
        self.alpha_slider.grid(row=6, column=1, pady=5)
        
        # System rotation control
        tk.Label(self.control_frame, text="System Rotation (°):").grid(row=7, column=0, sticky=tk.W, pady=3)
        self.phi_g_slider = tk.Scale(self.control_frame, from_=-180, to=180, orient=tk.HORIZONTAL, length=200,
                                    command=lambda x: self.update_parameter('phi_g', np.radians(float(x))))
        self.phi_g_slider.set(np.degrees(self.phi_g))
        self.phi_g_slider.grid(row=7, column=1, pady=5)
        
        # Angular velocity control
        tk.Label(self.control_frame, text="Angular Velocity (rad/s):").grid(row=8, column=0, sticky=tk.W, pady=3)
        self.omega_slider = tk.Scale(self.control_frame, from_=0.1, to=5.0, resolution=0.1, orient=tk.HORIZONTAL, length=200,
                                    command=lambda x: self.update_parameter('omega_1', float(x)))
        self.omega_slider.set(self.omega_1)
        self.omega_slider.grid(row=8, column=1, pady=5)
        
        # Configuration control
        self.config_var = tk.IntVar(value=1)
        tk.Radiobutton(self.control_frame, text="Configuration 1", variable=self.config_var, value=1,
                      command=lambda: self.update_parameter('config', 1)).grid(row=9, column=0, pady=5)
        tk.Radiobutton(self.control_frame, text="Configuration 2", variable=self.config_var, value=2,
                      command=lambda: self.update_parameter('config', 2)).grid(row=9, column=1, pady=5)

        
        self.play_button = tk.Button(self.control_frame, text="▶ Play", command=self.start_animation)
        self.play_button.grid(row=10, column=0, padx=5, pady=5)
        
        self.stop_button = tk.Button(self.control_frame, text="■ Stop", command=self.stop_animation)
        self.stop_button.grid(row=10, column=1, padx=5, pady=5)
        
        # Analysis mode selection
        analysis_frame = tk.LabelFrame(self.control_frame, text="Analysis View")
        analysis_frame.grid(row=11, column=0, columnspan=2, pady=10, padx=5, sticky=tk.W+tk.E)
        
        self.analysis_var = tk.StringVar(value="Mechanism")
        tk.Radiobutton(analysis_frame, text="Mechanism", variable=self.analysis_var, value="Mechanism",
                      command=lambda: self.set_analysis_mode("Mechanism")).grid(row=0, column=0, sticky=tk.W, pady=2)
        tk.Radiobutton(analysis_frame, text="Position Analysis", variable=self.analysis_var, value="Position",
                      command=lambda: self.set_analysis_mode("Position")).grid(row=0, column=1, sticky=tk.W, pady=2)
        tk.Radiobutton(analysis_frame, text="Velocity Analysis", variable=self.analysis_var, value="Velocity",
                      command=lambda: self.set_analysis_mode("Velocity")).grid(row=0, column=2, sticky=tk.W, pady=2)
        
        # Status display
        self.status_frame = tk.LabelFrame(self.control_frame, text="Mechanism Status")
        self.status_frame.grid(row=12, column=0, columnspan=2, pady=10, padx=5, sticky=tk.W+tk.E)
        
        self.grashof_label = tk.Label(self.status_frame, text="Checking Grashof condition...")
        self.grashof_label.grid(row=0, column=0, columnspan=2, pady=5)
        
        # Default mechanism presets
        preset_frame = tk.LabelFrame(self.control_frame, text="Preset Mechanisms")
        preset_frame.grid(row=13, column=0, columnspan=2, pady=10, padx=5, sticky=tk.W+tk.E)
        
        self.preset_var = tk.StringVar()
        self.preset_menu = ttk.Combobox(preset_frame, textvariable=self.preset_var, width=25)
        self.preset_menu['values'] = (
            "Grashof Double Crank",
            "Grashof Crank Rocker",
            "Grashof Double Rocker",
            "Non-Grashof 1",
            "Non-Grashof 2"
        )
        self.preset_menu.grid(row=0, column=0, padx=5, pady=5)
        self.preset_menu.current(0)
        
        self.load_preset_button = tk.Button(preset_frame, text="Load Preset", command=self.load_preset)
        self.load_preset_button.grid(row=0, column=1, padx=5, pady=5)

    def update_parameter(self, param, value):
        """Update a single parameter and refresh the display"""
        setattr(self, param, value)
        try:
            # Update the linkage parameters
            params = {
                'L0': self.L0,
                'L1': self.L1,
                'L2': self.L2,
                'L3': self.L3,
                'Lap': self.Lap,
                'alpha': self.alpha,
                'phi_g': self.phi_g,
                'config': self.config,
                'omega_1': self.omega_1
            }
            self.linkage.update_parameters(**params)
            
            # Update Grashof status
            grashof_status = self.linkage.Grashof()
            if grashof_status == "Grashof":
                self.grashof_label.config(text="Grashof Mechanism", fg="green")
            elif grashof_status == "Non-Grashof":
                self.grashof_label.config(text="Non-Grashof Mechanism", fg="orange")
            else:
                self.grashof_label.config(text=f"Invalid Mechanism", fg="red")
                self.animation_running = False
                # Reset the display or show error state
                self.fig.clear()
                ax = self.fig.add_subplot(111)
                ax.text(0.5, 0.5, 'Invalid Configuration\nPlease adjust parameters', 
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=ax.transAxes)
                self.canvas.draw()
                
            # Update display based on current analysis mode
            self.update_display()
            
        except Exception as e:
            print(f"Error updating linkage: {e}")
            self.grashof_label.config(text=f"Invalid Configuration: {str(e)}", fg="red")
            # Reset the display or show error state
            self.animation_running = False
            self.fig.clear()
            ax = self.fig.add_subplot(111)
            ax.text(0.5, 0.5, 'Invalid Configuration\nPlease adjust parameters', 
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            self.canvas.draw()

    def update_linkage(self):
        """Initialize or reinitialize the linkage"""
        try:
            self.linkage = FourBarLinkage()
            self.linkage.omega_1 = self.omega_1
            self.update_parameter('config', self.config)  # This will trigger a full update
        except Exception as e:
            print(f"Error creating linkage: {e}")
            self.grashof_label.config(text=f"Error: {str(e)}", fg="red")
    
    def set_analysis_mode(self, mode):
        """Change the analysis mode and update the display"""
        self.analysis_mode = mode
        self.update_display()
    
    def start_animation(self):
        """Start the animation"""
        if not self.animation_running and hasattr(self, 'linkage') and self.linkage.positions:
            self.animation_running = True
            self.animate(self.current_frame)
    
    def stop_animation(self):
        """Stop the animation"""
        self.animation_running = False
    
    def animate(self, frame_index):
        """Animate the mechanism"""
        if not self.animation_running:
            return
        
        if not hasattr(self, 'linkage') or not self.linkage.positions:
            self.animation_running = False
            return
            
        self.current_frame = frame_index
        
        if self.analysis_mode == "Mechanism":
            self.update_mechanism_display(frame_index)
        
        next_frame = (frame_index + 1) % len(self.linkage.positions)
        self.root.after(int(100 - self.animation_speed), lambda: self.animate(next_frame))
    
    def update_display(self):
        """Update the display based on current analysis mode"""
        if not hasattr(self, 'linkage') or not self.linkage.positions:
            return
            
        if self.analysis_mode == "Mechanism":
            self.update_mechanism_display(self.current_frame)
        elif self.analysis_mode == "Position":
            self.update_position_display()
        elif self.analysis_mode == "Velocity":
            self.update_velocity_display()
    
    def update_mechanism_display(self, frame_index):
        """Update the mechanism display"""
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        self.linkage.plot_simulation(ax=ax, index=frame_index)
        ax.grid(True)
        ax.set_aspect('equal')
        
        # Calculate the maximum extent of all trajectories for consistent scaling
        all_points = np.array(self.linkage.traj_A + self.linkage.traj_B + self.linkage.traj_P)
        if len(all_points) > 0:
            ax.set_ylabel('Y Position')
            ax.set_xlabel('X Position')
            ax.set_title("Four-Bar Linkage Simulation")

        
        self.canvas.draw()
    
    def update_position_display(self):
        """Update the position analysis display"""
        self.fig.clear()
        self.linkage.plot_position_graphs(self.fig)
        self.canvas.draw()
    
    def update_velocity_display(self):
        """Update the velocity analysis display"""
        self.fig.clear()
        self.linkage.plot_velocity_graphs(self.fig)
        self.canvas.draw()
    
    def load_preset(self):
        """Load a preset mechanism configuration"""
        preset = self.preset_var.get()
        
        presets = {
            # Grashof & ground L0 = shortest → Double‑Crank
            "Grashof Double Crank": {
                "L0": 12, "L1": 56, "L2": 67, "L3": 37,
                "Lap": 15, "alpha": np.radians(30), "phi_g": np.radians(0),
                "config": 1
            },

            # Grashof & ground L0 = intermediate → Crank‑Rocker
            "Grashof Crank Rocker": {
                "L0": 40, "L1": 10, "L2": 30, "L3": 50,
                # sorted: [10,30,40,50] → s+l = 10+50 = 60 < p+q = 30+40 = 70
                "Lap": 20, "alpha": np.radians(45), "phi_g": np.radians(45),
                "config": 1
            },

            # Grashof & ground L0 = longest → Double‑Rocker
            "Grashof Double Rocker": {
                "L0": 60, "L1": 40, "L2": 50, "L3": 20,
                # sorted: [20,50,55,60] → s+l = 20+60 = 60 < p+q = 50+55 = 105
                "Lap": 10, "alpha": np.radians(45), "phi_g": np.radians(45),
                "config": 1
            },

            # Non‑Grashof (s+l > p+q)
            "Non-Grashof 1": {
                "L0": 55, "L1": 45, "L2": 35, "L3": 30,
                # sorted: [30,35,45,55] → s+l = 30+55 = 85 > p+q = 35+45 = 80
                "Lap": 15, "alpha": np.radians(45), "phi_g": np.radians(45),
                "config": 1
            },

            # Non‑Grashof (already valid)
            "Non-Grashof 2": {
                "L0": 70, "L1": 50, "L2": 40, "L3": 93,
                # sorted: [40,50,70,93] → s+l = 40+93 = 133 > p+q = 50+70 = 120
                "Lap": 12, "alpha": np.radians(45), "phi_g": np.radians(45),
                "config": 1
            }
        }

        if preset in presets:
            config = presets[preset]
            
            # Update parameters
            self.L0 = config["L0"]
            self.L1 = config["L1"]
            self.L2 = config["L2"]
            self.L3 = config["L3"]
            self.Lap = config["Lap"]
            self.alpha = config["alpha"]
            self.phi_g = config["phi_g"]
            self.config = config["config"]
            
            # Update sliders
            self.L0_slider.set(self.L0)
            self.L1_slider.set(self.L1)
            self.L2_slider.set(self.L2)
            self.L3_slider.set(self.L3)
            self.Lap_slider.set(self.Lap)
            self.alpha_slider.set(np.degrees(self.alpha))
            self.phi_g_slider.set(np.degrees(self.phi_g))
            self.config_var.set(self.config)
            
            # Update the linkage
            self.update_linkage()


def show_simulation():
    root = tk.Tk()
    root.title("Four-Bar Linkage Simulation with Position and Velocity Analysis")
    simulator = FourBarLinkageSimulator(root)
    root.mainloop()

def show_license_info():
    license_text = (
        "MIT License\n\n"
        "Copyright (c) 2025 Amelia Hoyos\n\n"
        "Permission is hereby granted, free of charge, to any person obtaining a copy "
        "of this software and associated documentation files (the \"Software\"), to deal "
        "in the Software without restriction, including without limitation the rights "
        "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell "
        "copies of the Software, and to permit persons to whom the Software is "
        "furnished to do so, subject to the following conditions:\n\n"
        "The above copyright notice and this permission notice shall be included in "
        "all copies or substantial portions of the Software.\n\n"
        "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR "
        "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, "
        "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE "
        "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER "
        "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, "
        "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN "
        "THE SOFTWARE.\n\n"
        "Author: Amelia Hoyos"
    )
    messagebox.showinfo("License and Author", license_text)


if __name__ == "__main__":
    show_simulation()