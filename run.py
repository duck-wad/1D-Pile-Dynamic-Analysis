import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# File paths
displacement_file_path = "data\\OUTPUT_DISPLACEMENT.csv"
input_file_path = "INPUT.txt"
cpp_executable = "program.exe"  # Path to the compiled C++ executable

# Function to run the C++ program
def run_cpp_program(input_file_path):
    """Runs the C++ executable with the input file."""
    print("Running C++ program...")
    try:
        subprocess.run([cpp_executable, input_file_path], check=True)
        print("C++ program completed.")
    except FileNotFoundError:
        print(f"Error: C++ executable '{cpp_executable}' not found.")
        exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error: C++ program failed with error code {e.returncode}.")
        exit(1)

# Function to read pile length from INPUT.txt
def read_input_file(file_path):
    """Reads the pile length from the input file."""
    pile_length = 0
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("length:"):
                pile_length = float(line.split(":")[1].strip())
    return pile_length

# Run the C++ program before continuing
run_cpp_program(input_file_path)

# Read pile length
pile_length = read_input_file(input_file_path)

# Load displacement data from the CSV file
data = np.loadtxt(displacement_file_path, delimiter=",")
time = data[:, 0]  # Time values (first column)
displacement_data = data[:, 1:]  # Displacement values (remaining columns)

num_nodes = displacement_data.shape[1]  # Number of nodes
time_steps = displacement_data.shape[0]  # Number of time steps

# Determine the maximum displacement for scaling
max_displacement = np.max(np.abs(displacement_data))  # Largest absolute displacement
x_limit = max_displacement * 1.2  # Add some padding for better visibility

# Node positions along the length of the pile (flipped: 0 at the top)
node_positions = np.linspace(0, pile_length, num_nodes)  # Top of the pile at 0 m

# Plot settings
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(-x_limit, x_limit)  # Dynamic horizontal limits based on max displacement
ax.set_ylim(pile_length, 0)  # Flip the vertical axis: 0 m at top, pile length at bottom
ax.set_xlabel("Displacement (m)")
ax.set_ylabel("Pile Length (m)")
ax.set_title("Pile Movement Animation")

# Add a timer as a text annotation
time_text = ax.text(0.95, 0.05, "", transform=ax.transAxes, ha="right", va="bottom", fontsize=12)

# Initialize line objects for nodes
lines = ax.plot([], [], 'o-', lw=2, markersize=5, label="Nodes")
line = lines[0]

# Initialize animation
def init():
    """Initialize the line objects."""
    line.set_data([], [])
    time_text.set_text("")
    return line, time_text

def animate(i):
    """Update the positions of nodes and the timer at each time step."""
    x = displacement_data[i, :]  # Displacement at time step 'i' for all nodes
    y = node_positions           # Vertical positions remain fixed
    line.set_data(x, y)          # Update line data
    time_text.set_text(f"Time: {time[i]:.3f} s")  # Update the timer with correct time
    return line, time_text

# Create the animation
# The interval is set to time step in milliseconds (0.01 seconds = 10 ms)
real_time_interval = time[1] * 1000  # Convert time step to milliseconds
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=time_steps, interval=real_time_interval, blit=True)

# Show animation
plt.legend()
plt.show()
