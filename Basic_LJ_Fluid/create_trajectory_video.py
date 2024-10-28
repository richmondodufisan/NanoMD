import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def parse_lammps_dump(traj_file):
    frames = []
    with open(traj_file, 'r') as file:
        while True:
            # Read through each frame in the dump file
            line = file.readline()
            if not line:
                break  # End of file
            
            if "ITEM: TIMESTEP" in line:
                timestep = int(file.readline().strip())
                file.readline()  # Skip ITEM: NUMBER OF ATOMS
                num_atoms = int(file.readline().strip())
                file.readline()  # Skip ITEM: BOX BOUNDS
                file.readline()
                file.readline()
                file.readline()
                file.readline()  # Skip ITEM: ATOMS id type x y z

                # Initialize lists to store positions and types
                positions = []
                types = []
                
                for _ in range(num_atoms):
                    atom_data = file.readline().split()
                    types.append(int(atom_data[1]))
                    positions.append([float(atom_data[2]), float(atom_data[3]), float(atom_data[4])])

                frames.append((positions, types))
    return frames

def create_video_two_types(traj_file, output_file, dpi=300):
    # Parse the LAMMPS dump file
    frames = parse_lammps_dump(traj_file)
    
    # Only keep every Nth frame
    frames = frames[::1]

    # Set up the figure and 3D axis for the animation
    fig = plt.figure()
    fig.patch.set_facecolor('black')  # Set figure background to black
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')  # Set plot background to black

    # Get global limits for consistent scaling
    all_positions = [pos for frame in frames for pos in frame[0]]
    x_min, x_max = min(p[0] for p in all_positions), max(p[0] for p in all_positions)
    y_min, y_max = min(p[1] for p in all_positions), max(p[1] for p in all_positions)
    z_min, z_max = min(p[2] for p in all_positions), max(p[2] for p in all_positions)

    def update(frame_index):
        ax.cla()  # Clear the axis
        positions, atom_types = frames[frame_index]
        
        # Separate positions based on atom types
        pos_type1 = [positions[i] for i in range(len(atom_types)) if atom_types[i] == 1]
        pos_type2 = [positions[i] for i in range(len(atom_types)) if atom_types[i] == 2]

        # Scatter plot for each atom type with specified colors
        if pos_type1:
            ax.scatter(*zip(*pos_type1), color='blue', s=5.0, label="Type 1")
        if pos_type2:
            ax.scatter(*zip(*pos_type2), color='red', s=5.0, label="Type 2")
        
        # Set consistent limits for x, y, z axes
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_zlim(z_min, z_max)
        ax.set_box_aspect([x_max - x_min, y_max - y_min, z_max - z_min])  # Aspect ratio
        ax.set_title(f"Frame {frame_index * 10 + 1}", color='white')  # Frame number with white text
        ax.tick_params(colors='white')  # White tick labels for contrast
        ax.legend(loc="upper left")

    # Create the animation
    anim = FuncAnimation(fig, update, frames=len(frames), interval=50)

    # Save the animation as a high-resolution video file
    anim.save(output_file, fps=20, dpi=dpi)

# Usage
create_video_two_types("dump.lammpstrj", "output_video.mp4", dpi=300)
