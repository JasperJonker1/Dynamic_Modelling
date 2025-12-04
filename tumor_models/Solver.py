import models
from matplotlib import pyplot as plt


class Solver:
    def __init__(self, n_steps, time, model, init_volume=1e-7, init_time=0, *modelvalues):
        """

        :param n_steps: number of steps to take in the function
        :param time: amount of time (arbitrary value)
        :param init_volume: V0 value for the Solver
        :param model: model from model.py to use for the Solver
        :param modelvalues: __init__ values for the given model
        """
        self.n_steps = n_steps
        self.dt = time / self.n_steps
        self.v_zero = init_volume
        self.t_zero = init_time
        self.model = model(*modelvalues)
        self.algorithm = self.model.dVdt

    def euler_function(self):
        """
        Solve using the given algorithm with the euler function
        :return: list of time values and list of volume values at the time values
        """
        t = self.t_zero
        V = self.v_zero

        ts = [t]
        Vs = [V]

        for i in range(self.n_steps):
            dy1 = self.algorithm(V, t) * self.dt
            t = t + self.dt
            V = V + dy1

            ts.append(t)
            Vs.append(V)

        return ts, Vs

    def heun_function(self):
        """
        Solve using the given algorithm with the heun function
        :return: list of time values and list of volume values at the time values
        """
        t = self.t_zero
        V = self.v_zero

        ts = [t]
        Vs = [V]

        for i in range(self.n_steps):
            # first update
            dydt1 = self.algorithm(V, t)
            t1 = t + self.dt
            V1 = V + dydt1 * self.dt
            # second update using Heun method
            dydt2 = self.algorithm(V1, t1)
            t = t + self.dt
            V = V + (dydt1 + dydt2) / 2.0 * self.dt

            ts.append(t)
            Vs.append(V)

        return ts, Vs

    def runge_kutta_function(self):
        """
        Solve using the given algorithm using the standard runge-kutta function
        :return: list of time values and list of volume values at the time values
        """
        t = self.t_zero
        V = self.v_zero

        ts = [t]
        Vs = [V]

        for i in range(self.n_steps):
            dydt1 = self.algorithm(V, t)
            t1 = t + 0.5 * self.dt
            V1 = V + 0.5 * dydt1 * self.dt
            # second update
            dydt2 = self.algorithm(V1, t1)
            t2 = t + 0.5 * self.dt
            V2 = V + 0.5 * dydt2 * self.dt
            # third update
            dydt3 = self.algorithm(V2, t2)
            t3 = t + self.dt
            V3 = V + dydt2 * self.dt
            # fourth update using base Runge-Kutta
            dydt4 = self.algorithm(V3, t3)
            t = t + self.dt
            V = V + (dydt1 + 2.0 * dydt2 + 2.0 * dydt3 + dydt4) / 6.0 * self.dt

            ts.append(t)
            Vs.append(V)

        return ts, Vs

    def __str__(self):
        return f"n_steps: {self.n_steps}, dt: {self.dt}, model: {self.model}, V0: {self.v_zero}, t0: {self.t_zero}"


def main():
    solver = Solver(1000, 10, models.VonBertalanffyModel, 1e-7, 0, 1, 1)
    time, volume = solver.runge_kutta_function()

    plt.plot(time, volume)
    plt.show()


if __name__ == "__main__":
    main()
