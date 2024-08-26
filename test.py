import matplotlib.pyplot as plt
import numpy as np

# Degrees array from 0 to 360
degrees = np.linspace(0, 360, 1001)

def func(degree):
    # Convert degree to radians for computation
    radian = np.deg2rad(degree)
    
    # Term 1: Larger amplitude (global maximum)
    k_theta1 = 5   # Larger force constant for global maximum
    n1 = 3
    phi1 = np.deg2rad(0)  # No phase shift
    
    # Term 2: Smaller amplitude (local maximum)
    k_theta2 = 2   # Smaller force constant for local maximum
    n2 = 1
    phi2 = np.deg2rad(180)  # Phase shift of 180 degrees
    
    # Combined potential energy function
    energy = (k_theta1 * (1 + np.cos(n1 * radian - phi1)) +
              k_theta2 * (1 + np.cos(n2 * radian - phi2)))
    
    return energy

# Calculate potential energy for each degree value
pot_energy = [func(degree) for degree in degrees]

# Plotting the potential energy vs degrees
plt.plot(degrees, pot_energy)
plt.xlabel('Dihedral Angle (degrees)')
plt.ylabel('Potential Energy')
plt.title('Torsional Potential Energy with Global and Local Maxima')
plt.grid(True)
plt.savefig('ptenergyfunc.png')
plt.show()
