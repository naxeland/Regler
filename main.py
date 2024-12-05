from scipy.optimize import root_scalar
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Define the system (numerator and denominator of the transfer function)
num = [16, 8]
den = [1, 16, 32, 32, 0]
# Target value in radians (-45 degrees)
target_value_final = -3*np.pi / 4

# Define the function based on the given equation
def func_final(x):
    return (1/4) * np.arctan(2*x) - np.pi/2 - np.arctan(x/2) - 2*np.arctan(x/4) - target_value_final



# Define the frequency range for the Bode plot


# Create a transfer function system
system = signal.TransferFunction(num, den)


# Generate the Bode plot
w, mag, phase = signal.bode(system)
# Plot magnitude

plt.axvline(x=1, color='r', linestyle='--')  # Add vertical line at the first asymptote

plt.xlim(0.1, 10)
plt.semilogx(w, 10**(mag/20))
plt.title('Bode plot')
plt.xlabel('Frequency [rad/s]')
plt.ylabel('Magnitude [dB]')
plt.yscale('log')
plt.ylim(0.01, 15)
plt.grid(True, which='both')
plt.figure()

# Plot phase

plt.semilogx(w, phase)
plt.xlabel('Frequency [rad/s]')
plt.ylabel('Phase [degrees]')
plt.xlim(0.1, 10)
plt.ylim(-200, 90)
plt.plot(func_final(w))
plt.plot(2.038324649571132, -135, 'ro')
plt.plot(2.816950008351725, 0, 'bo')
plt.grid(True, which='both',)


f_interp = interp1d(phase, w)

# Desired y-value
desired_y = -135

# Get the corresponding x-value(s)
x_val = f_interp(desired_y)

print(f"x-value for y={desired_y}: {x_val}")

plt.show()


# Target value in radians (-135 degrees)
target_value_final = -3*np.pi / 4

# Define the function based on the given equation
def func_final(x):
    return np.arctan(2*x) - np.pi/2 - np.arctan(x/2) - 2*np.arctan(x/4) - target_value_final

# Solve using Brent's method within a reasonable bracket
solution_final = root_scalar(func_final, bracket=[-10, 10], method='brentq')

# Check if the solution was found
if solution_final.converged:
    print(f"Solution for x: {solution_final.root}")
else:
    print("No solution found.")
    
print(func_final(2.816950008351725))



