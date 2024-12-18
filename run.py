import subprocess
import numpy as np
import matplotlib.pyplot as plt

# -*- coding: utf-8 -*-

# File paths
displacement_file_path = "data\\OUTPUT_DISPLACEMENT.csv"
velocity_file_path = "data\\OUTPUT_VELOCITY.csv"
acceleration_file_path = "data\\OUTPUT_ACCELERATION.csv"
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

# Load displacement, velocity, and acceleration data from CSV files
displacement_data = np.loadtxt(displacement_file_path, delimiter=",")
velocity_data = np.loadtxt(velocity_file_path, delimiter=",")
acceleration_data = np.loadtxt(acceleration_file_path, delimiter=",")

# Extract time, displacement, velocity, and acceleration
time = displacement_data[:, 0]  # Time values (first column)
displacement_data = displacement_data[:, 1:]  # Displacement values (remaining columns)
velocity_data = velocity_data[:, 1:]  # Velocity values (remaining columns)
acceleration_data = acceleration_data[:, 1:]  # Acceleration values (remaining columns)

num_nodes = displacement_data.shape[1]  # Number of nodes
time_steps = displacement_data.shape[0]  # Number of time steps

# Plot settings for the stacked graphs
fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# Plot displacement for each node
for node in range(num_nodes):
    axes[0].plot(time, displacement_data[:, node], label=f"Node {node + 1}")

axes[0].set_ylabel("Displacement (m)")
axes[0].set_title("Displacement vs Time for Each Node")
axes[0].legend()
axes[0].grid(True)

# Plot velocity for each node
for node in range(num_nodes):
    axes[1].plot(time, velocity_data[:, node], label=f"Node {node + 1}")

axes[1].set_ylabel("Velocity (m/s)")
axes[1].set_title("Velocity vs Time for Each Node")
axes[1].legend()
axes[1].grid(True)

# Plot acceleration for each node
for node in range(num_nodes):
    axes[2].plot(time, acceleration_data[:, node], label=f"Node {node + 1}")

axes[2].set_xlabel("Time (s)")
axes[2].set_ylabel("Acceleration (m/s^2)")
axes[2].set_title("Acceleration vs Time for Each Node")
axes[2].legend()
axes[2].grid(True)

# Adjust layout for better spacing
plt.tight_layout()

# Show the stacked plots
plt.show()
