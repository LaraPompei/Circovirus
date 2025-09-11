import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, CubicSpline

# Load antibody data
# Assumes CSV has: time in first column, antibody levels in second column
data = np.genfromtxt('./dados/viremiaPorcoInoculado.csv', delimiter=',', skip_header=1)
t = data[:, 0]
antibody = data[:, 1]

# Build interpolation functions
linear_interp = interp1d(t, antibody, kind='linear', fill_value='extrapolate')
cubic_interp = CubicSpline(t, antibody)

# Make a fine time grid
t_fine = np.linspace(t.min(), t.max(), 500)
antibody_linear = linear_interp(t_fine)
antibody_cubic = cubic_interp(t_fine)

# Plot
plt.figure(figsize=(10, 6))
plt.scatter(t, antibody, color='red', label='Data points')
plt.plot(t_fine, antibody_linear, '--', label='Linear interpolation')
plt.plot(t_fine, antibody_cubic, label='Cubic spline interpolation')

plt.xlabel("Time (days)")
plt.ylabel("Antibody² (units)")
plt.title("Antibody² Interpolation")
plt.legend()
plt.grid(True)
plt.show()

