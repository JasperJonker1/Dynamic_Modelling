from typing import List, Tuple
from . import models          # ⬅ relatieve import!
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
        """
        Numerieke ODE-oplosser voor tumorgroeimodellen.

        Parameters
        ----------
        n_steps : int
            Aantal tijdstappen.
        time : float
            Totale tijd (van t0 tot t0 + time).
        model : class
            Modelklasse uit models.py (bijv. models.VonBertalanffyModel).
        init_volume : float
            Beginvolume V0.
        init_time : float
            Begintijd t0.
        modelvalues :
            Parameters voor de model-__init__ (bijv. c, d, Vmax, ...).
        """
        self.n_steps = n_steps
        self.dt = time / self.n_steps
        self.v_zero = init_volume
        self.t_zero = init_time
        self.model = model(*modelvalues)
        self.algorithm = self.model.dVdt

    def euler_function(self) -> Tuple[List[float], List[float]]:
        """
        Los de ODE op met de expliciete Euler-methode.

        Returns
        -------
        ts : list
            Tijdstappen.
        Vs : list
            Volume op elk tijdstip.
        """
        t = self.t_zero
        V = self.v_zero

        ts = [t]
        Vs = [V]

        for _ in range(self.n_steps):
            dV = self.algorithm(V, t) * self.dt
            t = t + self.dt
            V = V + dV

            ts.append(t)
            Vs.append(V)

        return ts, Vs

    def heun_function(self) -> Tuple[List[float], List[float]]:
        """
        Los de ODE op met de Heun-methode (verbeterde Euler).
        """
        t = self.t_zero
        V = self.v_zero

        ts = [t]
        Vs = [V]

        for _ in range(self.n_steps):
            # eerste helling
            dVdt1 = self.algorithm(V, t)
            t1 = t + self.dt
            V1 = V + dVdt1 * self.dt

            # tweede helling
            dVdt2 = self.algorithm(V1, t1)

            # update met gemiddelde helling
            t = t + self.dt
            V = V + 0.5 * (dVdt1 + dVdt2) * self.dt

            ts.append(t)
            Vs.append(V)

        return ts, Vs

    def runge_kutta_function(self) -> Tuple[List[float], List[float]]:
        """
        Los de ODE op met de klassieke Runge-Kutta 4-methode.
        """
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

            t = t + self.dt
            V = V + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0 * self.dt

            ts.append(t)
            Vs.append(V)

        return ts, Vs

    def __str__(self) -> str:
        return (
            f"Solver(n_steps={self.n_steps}, dt={self.dt}, "
            f"model={self.model}, V0={self.v_zero}, t0={self.t_zero})"
        )


def main():
    solver = Solver(1000, 10.0, models.VonBertalanffyModel, 1e-7, 0.0, 1.0, 1.0)
    time, volume = solver.runge_kutta_function()

    plt.plot(time, volume)
    plt.xlabel("t")
    plt.ylabel("V(t)")
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()

