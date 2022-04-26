import numpy as np
import matplotlib.pyplot as plt
import math

step = 0.01
T = 848
e_dot = 1
t_max = 1
b = 0.25e-9
D = 30
mu = 45000
Q = 238000
rho_0 = 10e4
R = 8.314
Z = e_dot * math.exp(Q / R / T)
rho = rho_0
t_0 = 0
t = t_0
# time_values = []
rho_values = []
t_cr = t_max


a1 = 2.1
a2 = 176.0
a3 = 19.5
a4 = 0.000148
a5 = 151.0
a6 = 0.973
a7 = 5.77
a8 = 1.0
a9 = 0.0
a10 = 0.262
a11 = 0.0
a12 = 0.000605
a13 = 0.167

tau = 10e6 * mu * b ** 2
l = a1 * 1e-3 / Z ** a13
A1 = 1 / b / l
A2 = a2 * e_dot ** (-a9) * math.exp(-a3 * 1e3 / R / T)
A3 = a4 * 3e10 * tau / D * math.exp(-a5 * 1e3 / R / T)
rho_cr = -a11 * 1e13 + a12 * 1e13 * Z ** a10
is_critical = False


def drho_dt():
  global is_critical
  if not is_critical:
    if rho > rho_cr:
      t_cr = t
      is_critical = True
      t2 = t - t_cr
      round(t2, 2)
      step_number = round(t2 / step)
      X = rho_values[step_number]
    else:
      X = 0
  else:
    t2 = t - t_cr
    math.round(t2, 2)
    step_number = math.round(t2 / step)
    X = rho_values[step_number]
  result = A1 * e_dot - A2 * e_dot * rho - A3 * rho ** a8 * X
  return result
  
  
while True:
  rho = rho + step * drho_dt()
  t += step
  rho_values.append(rho)
  # time_values.append(t)
  if t > t_max:
    break
print(rho_values)
