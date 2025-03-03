import ase.io
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def create_video(traj_file, output_file, particle_color='blue', particle_radius=100, dpi=300):
    # Load the LAMMPS trajectory file using ASE with correct format
    atoms = ase.io.read(traj_file, format='lammps-dump-text', index=':')
    
    # Only keep every 10th frame
    atoms = atoms[::10]
    
    # Set up the figure and 3D axis for the animation
    fig = plt.figure()
    fig.patch.set_facecolor('black')  # Set figure background to black
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')  # Set plot background to black

    # Get global limits for consistent scaling
    all_positions = [atom.get_positions() for atom in atoms]
    x_min = min(pos[:, 0].min() for pos in all_positions)
    x_max = max(pos[:, 0].max() for pos in all_positions)
    y_min = min(pos[:, 1].min() for pos in all_positions)
    y_max = max(pos[:, 1].max() for pos in all_positions)
    z_min = min(pos[:, 2].min() for pos in all_positions)
    z_max = max(pos[:, 2].max() for pos in all_positions)

    def update(frame):
        ax.cla()  # Clear the axis
        positions = atoms[frame].get_positions()
        ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                   color=particle_color, s=particle_radius)
        
        # Set consistent limits for x, y, z axes
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_zlim(z_min, z_max)
        ax.set_box_aspect([x_max - x_min, y_max - y_min, z_max - z_min])  # Aspect ratio
        ax.set_title(f"Frame {frame * 10 + 1}", color='white')  # Frame number with white text
        ax.tick_params(colors='white')  # White tick labels for contrast

    # Create the animation
    anim = FuncAnimation(fig, update, frames=len(atoms), interval=50)

    # Save the animation as a high-resolution video file
    anim.save(output_file, fps=20, extra_args=['-vcodec', 'libx264'], dpi=dpi)

# Usage
create_video("dump.lammpstrj", "output_video.mp4", particle_color='blue', particle_radius=1.0, dpi=300)
