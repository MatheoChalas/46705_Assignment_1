import numpy as np
import matplotlib.pyplot as plt

#Impedance relay settings for 3 zones
Zt1 = 0.8*(500/4500)*(7 +70j)
Zt2 = 1.2*(500/4500)*(7 +70j)
Zt3 = (500/4500)*(7 +70j)+1.2*(500/4500)*(8 +80j)


# Function to plot circle
def plot_circle(center, point_on_circle, radius):
    # Calculate coordenates of circle
    theta = np.linspace(0, 2*np.pi, 100) #Generation of 100 point around 360°
    x = center.real + radius * np.cos(theta)
    y = center.imag + radius * np.sin(theta)
    

    plt.plot(x, y)

# Calculate center of circles
center1 = Zt1/2
center2 = Zt2/2
center3 = Zt3/2

# Calculate radius of circle
radio2 = abs(Zt2)/2
radio3 = abs(Zt3)/2

# Plot circles
plt.figure(figsize=(8, 8))
plot_circle(center1,Zt1, abs(Zt1)/2)
plot_circle(center2,Zt2, radio2)
plot_circle(center3,Zt3, radio3)

# Shading differences between circle 3 and 2
area_diff_2_3 = np.pi * radio3**2 - np.pi * radio2**2

if area_diff_2_3 > 0:
    theta_diff_2_3 = np.linspace(0, 2*np.pi, 100)
    x_diff_2_3 = center3.real + radio3 * np.cos(theta_diff_2_3)
    y_diff_2_3 = center3.imag + radio3 * np.sin(theta_diff_2_3)
    plt.fill(x_diff_2_3, y_diff_2_3, 'lightblue', alpha=0.5,label='Zone 3')


# Shading circle 2
theta2 = np.linspace(0, 2*np.pi, 100)
x2 = center2.real + radio2 * np.cos(theta2)
y2 = center2.imag + radio2 * np.sin(theta2)
plt.fill(x2, y2, 'lightgray', alpha=0.5,label='Zone 2')

# Shading circle 1
theta1 = np.linspace(0, 2*np.pi, 100)
x1 = center1.real + (abs(Zt1)/2) * np.cos(theta1)
y1 = center1.imag + (abs(Zt1)/2) * np.sin(theta1)
plt.fill(x1, y1, 'gray', alpha=0.5,label='Zone 1')

# Ploting the values of settings (as points) in the complex plane
plt.scatter(Zt1.real,Zt1.imag, color='blue', label='Impedance 1')
plt.scatter(Zt2.real,Zt2.imag, color='orange', label='Impedance 2')
plt.scatter(Zt3.real, Zt3.imag, color='green', label='Impedance 3')


# Emergency point
plt.scatter(16.652,5.804, color='red', label='Emergency Point')


plt.xlabel('R')
plt.ylabel("X")
plt.title('Mho relay - Carachteristics')
plt.gca().set_aspect('equal', adjustable='box')


plt.grid(True)
plt.legend()


plt.show()