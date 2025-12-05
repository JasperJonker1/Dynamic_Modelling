
from typing import List, Tuple
from . import models
from matplotlib import pyplot as plt


class Solver:
    def __init__(
        self,
        n_steps: int,
        time: float,
        model,
        init_volume: float = 1e-7,
        init_time: float = 0.0,
        *modelvalues,
    ):
        self.n_steps = n_steps
        self.dt = time / self.n_steps
        self.v_zero = init_volume
        self.t_zero = init_time
        self.model = model(*modelvalues)
        self.algorithm = self.model.dVdt

    def euler_function(self) -> Tuple[List[float], List[float]]:
        t = self.t_zero
        V = self.v_zero
        ts = [t]
        Vs = [V]
        for _ in range(self.n_steps):
            dV = self.algorithm(V, t) * self.dt
            t += self.dt
            V += dV
            ts.append(t)
            Vs.append(V)
        return ts, Vs

    def heun_function(self) -> Tuple[List[float], List[float]]:
        t = self.t_zero
        V = self.v_zero
        ts = [t]
        Vs = [V]
        for _ in range(self.n_steps):
            k1 = self.algorithm(V, t)
            t1 = t + self.dt
            V1 = V + k1 * self.dt
            k2 = self.algorithm(V1, t1)
            t += self.dt
            V += 0.5 * (k1 + k2) * self.dt
            ts.append(t)
            Vs.append(V)
        return ts, Vs

    def runge_kutta_function(self) -> Tuple[List[float], List[float]]:
        t = self.t_zero
        V = self.v_zero
        ts = [t]
        Vs = [V]
        for _ in range(self.n_steps):
            k1 = self.algorithm(V, t)
            t1 = t + 0.5 * self.dt
            V1 = V + 0.5 * k1 * self.dt

            k2 = self.algorithm(V1, t1)
            t2 = t + 0.5 * self.dt
            V2 = V + 0.5 * k2 * self.dt

            k3 = self.algorithm(V2, t2)
            t3 = t + self.dt
            V3 = V + k3 * self.dt

            k4 = self.algorithm(V3, t3)

            t += self.dt
            V += (k1 + 2*k2 + 2*k3 + k4) / 6.0 * self.dt

            ts.append(t)
            Vs.append(V)
        return ts, Vs