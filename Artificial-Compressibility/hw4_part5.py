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
u_min_100 = [-0.18225488212325147, -0.20428819469127907, -0.21135904318042442]
v_min_100 = [-0.2110297208460899, -0.24088305226738468, -0.25067029055974616]
v_max_100 = [0.1498415790967516, 0.17017184781038114, 0.1772188544677155]

u_min_100_actual = -0.21090
v_min_100_actual = -0.24533
v_max_100_actual = 0.17527

"""
Re = 400
"""
u_min_400 = [-0.18262518113435913, -0.28083783153153635, -0.31506607492857447]
v_min_400 = [-0.27032763586139846, -0.3774396149842766, -0.43367168749250673]
v_max_400 = [0.16970628047725245, 0.2509997262917907, 0.2893171011770739]

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








