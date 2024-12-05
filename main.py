from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

# Define the system (numerator and denominator of the transfer function)
num = [16, 8]
den = [1, 16, 32, 32, 0]



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
plt.grid(True, which='both',)
plt.figure()

plt.show()