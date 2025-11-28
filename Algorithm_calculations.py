import math
from math import log
from matplotlib import pyplot as plt
from random import gauss

def von_bertalanffy(V, t):
    d = 1
    c = 1
    return (c * pow(V, 2 / 3)) - (d * V)


def gompertz_les(V, t):
    cap = 1000
    c = 1
    return c * V * log(cap / V)


def gompertz_paper(V, t):
    alpha = 1
    beta = 1
    return alpha * pow(math.e, (-beta * t)) * V


def euler_function(n_steps, time, init_volume, algorithm):
    t = 0
    V = init_volume

    ts = [t]
    Vs = [V]

    dt = time / n_steps

    for i in range(n_steps):
        dy1 = algorithm(V, t) * dt
        t = t + dt
        V = V + dy1

        ts.append(t)
        Vs.append(V)

    return ts, Vs


def heun_function(n_steps, time, init_volume, algorithm):
    t = 0
    V = init_volume

    ts = [t]
    Vs = [V]

    dt = time / n_steps

    for i in range(n_steps):
        # Eerste update
        dydt1 = algorithm(V, t)
        t1 = t + dt
        V1 = V + dydt1 * dt
        # Tweede update volgens Heun
        dydt2 = algorithm(V1, t1)
        t = t + dt
        V = V + (dydt1 + dydt2) / 2.0 * dt

        ts.append(t)
        Vs.append(V)

    return ts, Vs


def runge_kutta_function(n_steps, time, init_volume, algorithm):
    t = 0
    V = init_volume

    ts = [t]
    Vs = [V]

    dt = time / n_steps

    for i in range(n_steps):
        dydt1 = algorithm(V, t)
        t1 = t + 0.5 * dt
        V1 = V + 0.5 * dydt1 * dt
        # Tweede update
        dydt2 = algorithm(V1, t1)
        t2 = t + 0.5 * dt
        V2 = V + 0.5 * dydt2 * dt
        # Derde update
        dydt3 = algorithm(V2, t2)
        t3 = t + dt
        V3 = V + dydt2 * dt
        # Vierde update volgens Runge-Kutta
        dydt4 = algorithm(V3, t3)
        t = t + dt
        V = V + (dydt1 + 2.0 * dydt2 + 2.0 * dydt3 + dydt4) / 6.0 * dt

        ts.append(t)
        Vs.append(V)

    return ts, Vs


def mean_squared_error(ts, V_real, V_predict):
    # ts, V_predict = euler_function(1000, 15, 1e-7, gompertz_paper)
    mse = sum((predict_val - real_val) ** 2 for predict_val, real_val in zip(V_predict, V_real))
    mse = mse / (len(ts))
    return mse


def random_search(n_params, V_real):
    params = {i: 0.0 for i in range(n_params)}

    ts, Vs = runge_kutta_function(1000, 15, 1e-7, gompertz_paper)
    error_rate = mean_squared_error(ts, V_real=V_real, V_predict=Vs)

    for i in range(n_params):
        tries = 0
        while True:
            new_param = params[i] + gauss(mu=0.0, sigma=1.0)
            new_error_rate = mean_squared_error(ts, V_real=V_real, V_predict=Vs)
            if new_error_rate < error_rate:
                error_rate = new_error_rate
                tries = 0
            if error_rate < 1e-3 or tries >= 1000:
                break


def pattern_search(n_params, V_real):
    params = {i: 0.0 for i in range(n_params)}

    ts, Vs = runge_kutta_function(1000, 15, 1e-7, gompertz_paper)
    error_rate = mean_squared_error(ts, V_real=V_real, V_predict=Vs)

    for i in range(n_params):
        step_size = 5.0
        tries = 0
        while True:
            new_param = params[i] + step_size
            new_error_rate = mean_squared_error(ts, V_real=V_real, V_predict=Vs)
            if new_error_rate < error_rate:
                error_rate = new_error_rate
                tries = 0
            else:
                step_size *= -0.9
                tries += 1
            if error_rate < 1e-3 or tries >= 10000:
                break



def main():
    time, volume = runge_kutta_function(1000, 15, 1e-7, gompertz_paper)

    plt.plot(time, volume)
    plt.show()


if __name__ == "__main__":
    main()
