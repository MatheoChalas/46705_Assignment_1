import numpy as np
import matplotlib.pyplot as plt

# Impedance relay settings for 3 zones
Zt1 = 0.8*(500/4500)*(7 +70j)
Zt2 = 1.2*(500/4500)*(7 +70j)
Zt3 = (500/4500)*(7 +70j)+1.2*(500/4500)*(8 +80j)

# Functions to plot semicircles with rotation
def plot_rotated_semicircles(Z_values):
    for Z in Z_values:
        r = np.abs(Z)
        angle_Z = np.angle(Zt3)  # Zt3 angle
        rotation_angle = -(np.pi/2 - angle_Z)  # Rotation angle
        theta = np.linspace(0, np.pi, 100)
        x = r * np.cos(theta + rotation_angle).real
        y = r * np.sin(theta + rotation_angle).real
        plt.plot(x, y)

# Plot semcircunferences with rotation
plot_rotated_semicircles([Zt1, Zt2, Zt3])

# Calculate the angle of the line
angle_line = -(np.pi/2 - np.angle(Zt3))

# Coordinates of line
x_line = np.linspace(-20, 20, 100)
y_line = np.tan(angle_line) * x_line

# Plot line
plt.plot(x_line, y_line, 'k--')

# Ploting the values of settings (as points) in the complex plane
plt.scatter(Zt1.real,Zt1.imag, color='blue', label='Impedance 1')
plt.scatter(Zt2.real,Zt2.imag, color='orange', label='Impedance 2')
plt.scatter(Zt3.real, Zt3.imag, color='green', label='Impedance 3')

# Emergency point
plt.scatter(16.652,5.804, color='red', label='Emergency Point')


plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('R')
plt.ylabel('X')
plt.title('Relay with directional restraint')
plt.grid(True)


plt.xlim(-20, 30)
plt.ylim(-5, 30)

plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.legend()
plt.show()

