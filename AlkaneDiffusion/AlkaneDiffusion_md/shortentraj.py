import os

def shorten_xtc(alkane, alkane_atoms, alkane_size):
    """
    Shortens the trajectory file for a given molecule and size.
    """
    # Format the command string with the appropriate variables
    command_cut = f"gmx trjconv -f nvt2_{alkane}_{alkane_size}.xtc -b 0 -e 2000 -o nvt2_cut_{alkane}_{alkane_size}.xtc"
    os.system(command_cut)

if __name__ == '__main__':
    molecules = {17: 'pentane', 20: 'hexane', 23: 'heptane', 26: 'octane', 32: 'decane', 47: 'pentadecane'}
    sizes = [512, 1024, 2048]

    # Iterate over molecules and sizes
    for atoms, name in molecules.items():  # Corrected from .list() to .items()
        for size in sizes:
            shorten_xtc(name, atoms, size)
