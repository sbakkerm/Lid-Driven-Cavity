# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 11:56:17 2021

@author: asus
"""
import numpy as np
import matplotlib.pyplot as plt

#Re_vals = [100, 400]
N_vals = [16, 32, 64]

# u_min, v_min, v_max are ROWS
# 32, 64, 128 are columns
vals_100 = np.zeros((3,3))

x_arr = np.linspace(10,70,100)
y_arr = np.ones(100)

"""
Re = 100
"""
u_min_100 = [-0.19019791421148824, -0.20768335418703077, -0.21253553108865303]
v_min_100 = [-0.2271356601720124, -0.24766832604588787, -0.252445384744282]
v_max_100 = [0.15910391374593466, 0.17418819468999072, 0.17839730686823016]

u_min_100_actual = -0.21090
v_min_100_actual = -0.24533
v_max_100_actual = 0.17527

"""
Re = 400
"""
u_min_400 = [-0.22823468886280585, -0.29129275400054644, -0.31724510651526616]
v_min_400 = [-0.3322239409474378, -0.40403619647804623, -0.4383164666915708]
v_max_400 = [0.2151184137768686, 0.26567635453186117, 0.2906739695116219]

u_min_400_actual = -0.32726
v_min_400_actual = -0.44993
v_max_400_actual = 0.30203

"""
u_min plot
"""
plt.figure(dpi=800)
plt.plot(N_vals, u_min_100, 'g--.', label = "Re = 100")
plt.plot(x_arr, y_arr*u_min_100_actual, 'k')
plt.plot(N_vals, u_min_400, 'm--.', label = "Re = 400")
plt.plot(x_arr, y_arr*u_min_400_actual, 'k', label = "Actual")
plt.xlim(10, 70)
plt.title(r"$u_{min}$")
plt.legend()
plt.xticks(N_vals)
plt.xlabel("N")
plt.ylabel("Velocity at vertical centerline")
plt.show()

"""
v_min plot
"""
plt.figure(dpi=800)
plt.plot(N_vals, v_min_100, 'g--.', label = "Re = 100")
plt.plot(x_arr, y_arr*v_min_100_actual, 'k')
plt.plot(N_vals, v_min_400, 'm--.', label = "Re = 400")
plt.plot(x_arr, y_arr*v_min_400_actual, 'k', label = "Actual")
plt.xlim(10, 70)
plt.title(r"$v_{min}$")
plt.legend()
plt.xticks(N_vals)
plt.xlabel("N")
plt.ylabel("Velocity at horizontal centerline")
plt.show()

"""
v_max plot
"""
plt.figure(dpi=800)
plt.plot(N_vals, v_max_100, 'g--.', label = "Re = 100")
plt.plot(x_arr, y_arr*v_max_100_actual, 'k')
plt.plot(N_vals, v_max_400, 'm--.', label = "Re = 400")
plt.plot(x_arr, y_arr*v_max_400_actual, 'k', label = "Actual")
plt.xlim(10, 70)
plt.title(r"$v_{max}$")
plt.legend()
plt.xticks(N_vals)
plt.xlabel("N")
plt.ylabel("Velocity at horizontal centerline")
plt.show()








