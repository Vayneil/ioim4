import numpy as np
import matplotlib.pyplot as plt
import math

t_cr = None
is_critical = False
fig = plt.figure()
plt.xlabel("t, czas")
plt.ylabel("ro")


def xd(T=898.0, t_max=1.0, e_dot=1.0):
    global t_cr
    step = 0.001
    # e_dot = 1.0
    # t_max = 1.0
    # T = 898.0
    b = 0.25e-9
    D = 30.0
    mu = 45000.0
    Q = 238000.0
    rho_0 = 1e4
    R = 8.314
    Z = e_dot * math.exp(Q / R / T)
    rho = rho_0
    t_0 = 0.0
    t = t_0
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

    tau = 1e6 * mu * b ** 2 * 0.5
    l = a1 * 1e-3 / (Z ** a13)
    A1 = 1 / b / l
    A2 = a2 * e_dot ** (-a9) * math.exp(-a3 * 1e3 / R / T)
    A3 = a4 * 3e10 * tau / D * math.exp(-a5 * 1e3 / R / T)
    rho_cr = -a11 * 1e13 + a12 * 1e13 * Z ** a10
    is_critical = False

    def drho_dt():
        global is_critical
        global t_cr
        if not is_critical:
            if rho > rho_cr:
                t_cr = t
                is_critical = True
                t2 = t - t_cr
                round(t2, 3)
                step_number = round(t2 / step)
                X = rho_values[step_number]
            else:
                X = 0
        else:
            t2 = t - t_cr
            round(t2, 3)
            step_number = round(t2 / step)
            X = rho_values[step_number]
        result = A1 * e_dot - A2 * e_dot * rho - A3 * rho ** a8 * X
        return result


    while True:
        rho = rho + step * drho_dt()
        t += step
        rho_values.append(rho)
        if t > t_max:
            break
    # print(rho_values)
    plt.plot(np.arange(t_0, t_max, t_max / len(rho_values)), rho_values, label="T = " + str(T))
    # plt.show()

if __name__ == '__main__':
    for temp in [848.0, 873.0, 898.0]:
        xd(T=temp, t_max=0.1, e_dot=10)
        t_cr = None
        is_critical = False
    plt.legend()
    plt.show()
