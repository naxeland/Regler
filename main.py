from scipy.optimize import root_scalar
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from sympy import symbols, solve, Eq

# Define the system (numerator and denominator of the transfer function)
num = [0.5, 0.25]
den = [1/32, 1/2, 1, 1, 0]

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
plt.ylim(0.008, 14)
plt.grid(True, which='both')
plt.figure()

# Plot phase
plt.semilogx(w, phase)
plt.xlabel('Frequency [rad/s]')
plt.ylabel('Phase [degrees]')
plt.xlim(0.08, 12)
plt.ylim(-210, 170)
plt.plot(2.038324649571132, -135, 'ro')
plt.plot(2.816950008351725, -135, 'bo')
plt.grid(True, which='both',)

plt.show()


# Target value in radians (-135 degrees)
target_value_final = -3*np.pi / 4
# Target value in radians (-135 degrees)
target_value_final_pi =-np.pi

# Define the function based on the given equation
# Target value in radians (-135 degrees)
target_value_final = -3*np.pi / 4

# Define the function based on the given equation
def func_final(x):
    return (1/4) * np.arctan(2*x) - np.pi/2 - np.arctan(x/2) - 2*np.arctan(x/4) - target_value_final

def func_G(x):
    y = (1/4) * np.sqrt(4 * (x**2)) / (x * np.sqrt(1 + (x/2)**2) * np.sqrt(1 + (x/4)**2) * np.sqrt(1 + (x/4)**2))
    return y

def func_Kp(x):
    y = 1/func_G(x)
    return y

def func_wc(x):
    return np.arctan(2*x) - np.pi/2 - np.arctan(x/2) - 2*np.arctan(x/4) - target_value_final

def func_wpi(x):
    return np.arctan(2*x) - np.pi/2 - np.arctan(x/2) - 2*np.arctan(x/4) - target_value_final_pi

# Solve using Brent's method within a reasonable bracket
solution_final = root_scalar(func_wc, bracket=[-10, 10], method='brentq')

# Solve using Brent's method within a reasonable bracket
solution_final_pi = root_scalar(func_wpi, bracket=[-10, 10], method='brentq')

# Check if the solution was found
if solution_final.converged:
    print(f"Solution for Wc: {solution_final.root}")
else:
    print("No solution found.")
    
print('Solution for Wpi: ', solution_final_pi.root)
    
#print('Y value is: ', func_wc(2.816950008351725))

Kp = func_Kp(2.816950008351725)

print('Since Kp = 1/|G(x)|, then Kp: ', Kp)

wpi = func_G(solution_final_pi.root)

print('Using ', solution_final_pi.root ,'in |G(x)| results in: ', wpi)

L = func_G(solution_final_pi.root) * func_Kp(2.816950008351725)

print('L = |G(x)| * Kp, then L: ', L)

Am = 1/(func_G(solution_final_pi.root) * func_Kp(2.816950008351725))

print('Am = 1 / L, then Am: ', Am)

Amdb = 20 * np.log10(Am)

print('Amdb = 20log(Am), then Amdb: ', Amdb)



