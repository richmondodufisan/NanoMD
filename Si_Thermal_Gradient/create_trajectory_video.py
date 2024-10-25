import ase.io
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def create_video(traj_file, output_file, particle_color='blue', particle_radius=100):
    # Load the LAMMPS trajectory file using ASE with correct format
    atoms = ase.io.read(traj_file, format='lammps-dump-text', index=':')
    
    # Only keep every 100th frame
    atoms = atoms[::100]
    
    # Set up the figure and 3D axis for the animation
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio for x, y, z

    # Define function to update the plot for each frame
    def update(frame):
        ax.cla()  # Clear the axis
        positions = atoms[frame].get_positions()
        ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                   color=particle_color, s=particle_radius)
        ax.set_xlim([positions[:, 0].min(), positions[:, 0].max()])
        ax.set_ylim([positions[:, 1].min(), positions[:, 1].max()])
        ax.set_zlim([positions[:, 2].min(), positions[:, 2].max()])
        ax.set_title(f"Frame {frame * 100 + 1}")

    # Create the animation
    anim = FuncAnimation(fig, update, frames=len(atoms), interval=50)

    # Save the animation as a video file
    anim.save(output_file, fps=20, extra_args=['-vcodec', 'libx264'])

# Usage
create_video("dump.lammpstrj", "output_video.mp4", particle_color='red', particle_radius=1.0)
