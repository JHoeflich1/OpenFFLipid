import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt


u = mda.Universe('water_512.gro', 'nvt2_short_water_512_nojump.xtc')


atoms = u.select_atoms('all')


n_frames = len(u.trajectory)
n_atoms = len(atoms)
positions = np.zeros((n_frames, n_atoms, 3))  # Store the position of each atom in each frame


for i, ts in enumerate(u.trajectory):
    positions[i] = atoms.positions


def compute_msd(positions):
    """Calculate the mean squared displacement (MSD) for each atom."""
    n_frames = positions.shape[0]
    n_atoms = positions.shape[1]
    
    msd_per_atom = np.zeros((n_frames, n_atoms))
    
    for i in range(n_frames):
        # Calculate the displacement from the initial position
        displacement = positions[i:] - positions[:n_frames-i]
        squared_displacement = np.square(displacement).sum(axis=2)  # Sum over x, y, z
        msd_per_atom[i] = squared_displacement.mean(axis=0)  # Mean over frames for each atom

    return msd_per_atom

# Calculate MSD for each atom
msd_per_atom = compute_msd(positions)

# Convert MSD from Å² to nm²
msd_per_atom *= 0.01

# Calculate average MSD across all atoms
average_msd = msd_per_atom.mean(axis=1)

# Plotting MSDs
plt.figure(figsize=(10, 6))

# Plot MSD for each atom
for atom in range(n_atoms):
    plt.plot(msd_per_atom[:, atom], alpha=0.5, lw=0.5, label=f'Atom {atom+1}')

# Plot average MSD
plt.plot(average_msd, color='black', lw=2, label='Average MSD', linestyle='--')

plt.xlabel('Time Frame')
plt.ylabel('Mean Squared Displacement (nm²)')
plt.title('MSD for Each Atom and Average MSD')
plt.tight_layout()
plt.savefig('manualMSD.png')
