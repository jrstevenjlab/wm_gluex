# since mathematica and c++ aren't really suited for creating the
# ambiguity plots, I'll try using python, though its been forever
# since using this language, so it won't be pretty

from math import pi,cos, sqrt, atan
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# begin by asking user for input of 4 waves
print("\nEnter the magnitude^2 of your waves and phase "
      "differences (in radians).\n"
      "Parameters can be entered as a formula i.e."
      "3*pi/2\n"
      "These can be from the fit parameters output,\n"
      "or a set of generated parameters that you'd like to "
      "test for ambiguity\n")
"""
rho_pp = float(input("(rho_pp)^2: "))
rho_pm = float(input("(rho_pm)^2: "))
rho_mp = float(input("(rho_mp)^2: "))
rho_mm = float(input("(rho_mm)^2: "))
eta_p  = eval(input("eta_p: "))
eta_m  = eval(input("eta_m: "))
"""
rho_pp = 400.0
rho_pm = 225.0
rho_mp = 100.0
rho_mm = 169.0
eta_p  = atan(8/6)
eta_m  = atan(12/5)

K 	= rho_pp + rho_pm + rho_mp + rho_mm
beta_p  = rho_pp - rho_pm
beta_m  = rho_mp - rho_mm
Y_Phi   = 2*(sqrt(rho_pp*rho_mp)*cos(eta_p)
	  - sqrt(rho_pm*rho_mm)*cos(eta_m))
Y_phi   = 2*(sqrt(rho_pp*rho_mp)*cos(eta_p) 
	  + sqrt(rho_pm*rho_mm)*cos(eta_m))

print("Moments:\n"
      "K: \t%.2d\n"
      "B1: \t%.2d\nB-1: \t%.2f\n"
      "Y_Phi: \t%.2f\nY_phi: \t%.2f\n"
      % (K, beta_p, beta_m, Y_Phi, Y_phi))

rho_pp_min = max(beta_p, 0)
rho_pp_max = min(max((K+beta_m+beta_p)/2, 0),
                 max((K-beta_m+beta_p)/2, 0))

# put case here if inputs are non-ambiguous

# find rho bounds from bounds of rho_pp
rho_pm_min = max(rho_pp_min - beta_p, 0)
rho_pm_max = max(rho_pp_max - beta_p, 0)
rho_mp_min = max((K+beta_m+beta_p)/2 - rho_pp_max, 0)
rho_mp_max = max((K+beta_m+beta_p)/2 - rho_pp_min, 0)
rho_mm_min = max((K-beta_m+beta_p)/2 - rho_pp_max, 0)
rho_mm_max = max((K-beta_m+beta_p)/2 - rho_pp_min, 0)

print("Boundaries of rho_[m][e] values are:\n"
      "\t %.2f < (rho_pp)^2 < %.2f\n"
      "\t %.2f < (rho_pm)^2 < %.2f\n"
      "\t %.2f < (rho_mp)^2 < %.2f\n"
      "\t %.2f < (rho_mm)^2 < %.2f\n"
      % (rho_pp_min, rho_pp_max,
         rho_pm_min, rho_pm_max,
         rho_mp_min, rho_mp_max,
         rho_mm_min, rho_mm_max))

print("Normalized Boundaries:\n"
      "\t %.2f < (rho_pp)^2 < %.2f\n"
      "\t %.2f < (rho_pm)^2 < %.2f\n"
      "\t %.2f < (rho_mp)^2 < %.2f\n"
      "\t %.2f < (rho_mm)^2 < %.2f\n"
      % (rho_pp_min/K, rho_pp_max/K,
         rho_pm_min/K, rho_pm_max/K,
         rho_mp_min/K, rho_mp_max/K,
         rho_mm_min/K, rho_mm_max/K))


def phase_p(x):
    return (((Y_phi+Y_Phi)/4) / 
            (sqrt(x)*sqrt((K+beta_m+beta_p)/2-x)))

def phase_m(x):
    return (((Y_phi-Y_Phi)/4) / 
            (sqrt(x-beta_p)*sqrt((K-beta_m+beta_p)/2-x)))

# scan the cos(eta_p) and cos(eta_m) values for different rho_pp's
z = np.linspace(rho_pp_min, rho_pp_max)
x = [phase_p(i) for i in z[1:-1]]
y = [phase_m(i) for i in z[1:-1]]

# at bounds where phase is undefined, just define a line there
x_und = np.arange(-1.0, 1.5, 0.5)
y_und = np.arange(-1.0, 1.5, 0.5)

# handle our 3 different cases
case = [{"x":x, "y":y, "z":z[1:-1]}, 
        {"x":np.repeat(phase_p(z[0]), y_und.size), 
         "y":y_und, 
         "z":np.repeat(z[0], y_und.size)}, 
        {"x":x_und, 
         "y":np.repeat(phase_m(z[-1]), x_und.size), 
         "z":np.repeat(z[-1], x_und.size)}] 

# Make the plot
plt.rc('axes', labelsize=18)
plt.figure(1)
ax = plt.axes(projection='3d')

for dataset in case:
    ax.plot3D(dataset["x"], dataset["y"], dataset["z"], "gray")
#ax.plot3D(x, y, z[1:-1], 'gray')
ax.set_xlabel(r'$\cos(\eta^+)$')
ax.set_ylabel(r'$\cos(\eta^-)$')
ax.set_zlabel(r'$(\rho_1^+)^2$')


ax.set_xticks(np.arange(-1.0, 1.5, 0.5))
ax.set_yticks(np.arange(-1.0, 1.5, 0.5))

ax.set_xlim([1,-1])
ax.set_ylim([-1,1])

plt.savefig('ambiguityLine.pdf')

# Plot 2d Projection of phases
plt.rc('axes', labelsize=24)
plt.figure(2, figsize=(10,7))
ax_projection = plt.axes()
ax_projection.plot(case[0]["x"], case[0]["y"], color="gray")
ax_projection.plot(case[1]["x"], case[1]["y"], color="gray", 
                   linestyle="dashed")
ax_projection.plot(case[2]["x"], case[2]["y"], color="gray", 
                   linestyle="dashed",
                   label="Boundary Case")

ax_projection.set_xlim([0.4,1])
ax_projection.set_ylim([0.35,1])
ax_projection.set_xlabel(r'$\cos(\eta^+)$')
ax_projection.set_ylabel(r'$\cos(\eta^-)$')
ax_projection.xaxis.set_label_coords(0.89,-0.07)
ax_projection.yaxis.set_label_coords(-0.1,0.87)

ax_projection.grid()
minor_xticks = np.arange(0.4,1.02, 0.02)
minor_yticks = np.arange(0.36, 1.02, 0.02)
ax_projection.tick_params(axis="both", which="major", labelsize=16)
ax_projection.set_xticks(minor_xticks, minor=True)
ax_projection.set_yticks(minor_yticks, minor=True)

ax_projection.legend(fontsize="x-large")

plt.savefig("ambiguityProjection.pdf")
