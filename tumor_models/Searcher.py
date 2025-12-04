from random import gauss
from typing import Sequence, Dict, Optional
from .Solver import Solver
from . import models          
from matplotlib import pyplot as plt
import numpy as np


class Searcher:
    def __init__(
        self,
        real_vals: Sequence[float],
        time_values: Sequence[float],
        n_params: int,
        model,
        predict_function: Optional[str] = None,
    ):
        """
        Parameters
        ----------
        real_vals : lijst van gemeten V(t)
        time_values : lijst van tijdstippen
        n_params : aantal te optimaliseren parameters
        model : modelklasse (bijv. models.VonBertalanffyModel)
        predict_function : 'runge-kutta', 'heun' of None (Euler)
        """
        self.real_values = list(real_vals)
        self.time_values = list(time_values)
        self.n_params = n_params
        self.model = model
        self.predictor = predict_function

    def mean_squared_error(self, *params):
        

        # 1. Solver run met voldoende resolutie
        try:
            solver = Solver(
                2000,  # meer stappen → nauwkeurige tijd-as
                (self.time_values[-1] - self.time_values[0]),
                self.model,
                self.real_values[0],
                self.time_values[0],
                *params
            )

            # Kies methode
            if self.predictor == "runge-kutta":
                t_pred, V_pred = solver.runge_kutta_function()
            elif self.predictor == "heun":
                t_pred, V_pred = solver.heun_function()
            else:
                t_pred, V_pred = solver.euler_function()

        except:
            return float("inf")

        # 2. INTERPOLATIE naar de echte tijdspunten
        try:
            V_interp = np.interp(self.time_values, t_pred, V_pred)
        except:
            return float("inf")

        # 3. Log-MSE berekenen
        sumsq = 0.0
        for p, r in zip(V_interp, self.real_values):
            if p <= 0: p = 1e-9
            if r <= 0: r = 1e-9
            diff = np.log(p) - np.log(r)
            sumsq += diff * diff

        return sumsq / len(self.real_values)


    
    # ------------------------------------------------------------
    # 1. Random search per parameter
    # ------------------------------------------------------------

    def random_search_one(self, sigma: float = 1.0, max_tries: int = 10_000) -> Dict[int, float]:
        """
        Optimaliseer parameters één voor één door random updates.
        :return: dict {index: waarde}
        """
        params = {i: 0.0 for i in range(self.n_params)}
        error_rate = self.mean_squared_error(*params.values())

        for i in range(self.n_params):
            tries = 0
            while tries < max_tries:
                new_params = params.copy()
                # maak nieuwe parameter via een Gauss-stap
                new_params[i] = params[i] + gauss(0.0, sigma)
                new_error_rate = self.mean_squared_error(*new_params.values())

                if new_error_rate < error_rate:
                    params[i] = new_params[i]
                    error_rate = new_error_rate
                    tries = 0  # reset als we beter zijn
                else:
                    tries += 1
        return params

    # ------------------------------------------------------------
    # 2. Pattern search per parameter
    # ------------------------------------------------------------

    def pattern_search(self, initial_step: float = 10.0, tol: float = 1e-6) -> Dict[int, float]:
        """
        Eenvoudige pattern search per parameter.
        :return: dict met optimale parameters
        """
        params = {i: 0.0 for i in range(self.n_params)}
        error_rate = self.mean_squared_error(*params.values())

        for i in range(self.n_params):
            step_size = initial_step
            while abs(step_size) > tol:
                new_params = params.copy()
                new_params[i] = params[i] + step_size
                new_error_rate = self.mean_squared_error(*new_params.values())

                if new_error_rate < error_rate:
                    params[i] = new_params[i]
                    error_rate = new_error_rate
                else:
                    step_size *= -0.8  # keer om en maak kleiner
        return params

    # ------------------------------------------------------------
    # 3. Random search over alle parameters tegelijk
    # ------------------------------------------------------------

    def random_search_all(
        self,
        sigma: float = 0.01,
        max_tries: int = 1000,
    ) -> dict[int, float]:
        """
        Hill-climbing random search: alle parameters tegelijk updaten.

        :param sigma: standaardafwijking van de Gauss-stap.
        :param max_tries: max aantal opeenvolgende mislukte pogingen.
        :return: dict met gevonden parameters {index: waarde}
        """
        # start met alle parameters op 0.0
        params = {i: 0.0 for i in range(self.n_params)}
        error_rate = self.mean_squared_error(*params.values())

        # extra veiligheid: als hier al iets misgaat, meteen stoppen
        if error_rate is None or not isinstance(error_rate, (int, float)):
            error_rate = float("inf")

        tries = 0
        while tries < max_tries:
            new_params = {
                key: val + gauss(0.0, sigma)  # kleine gauss-stap rond huidige waarde
                for key, val in params.items()
            }
            new_error_rate = self.mean_squared_error(*new_params.values())

            # als er iets misgaat in MSE → beschouw als heel slecht
            if new_error_rate is None or not isinstance(new_error_rate, (int, float)):
                new_error_rate = float("inf")

            if new_error_rate < error_rate:
                params = new_params
                error_rate = new_error_rate
                tries = 0  # reset bij verbetering
            else:
                tries += 1

        return params

    # ------------------------------------------------------------
    # 4. Pattern search over alle parameters tegelijk
    # ------------------------------------------------------------

    def pattern_search_all(
        self,
        initial_step: float = 10.0,
        tol: float = 1e-6,
    ) -> Dict[int, float]:
        """
        Pattern search over alle parameters tegelijk.

        :return: dict met optimale parameters
        """
        params = {i: 0.0 for i in range(self.n_params)}
        step_size = {i: initial_step for i in range(self.n_params)}
        error_rate = self.mean_squared_error(*params.values())

        # Blijf zoeken zolang er minstens één stap nog groot is
        while max(abs(s) for s in step_size.values()) > tol:
            for key in params:
                new_params = params.copy()
                # eerst vooruit
                new_params[key] = params[key] + step_size[key]
                new_error_rate = self.mean_squared_error(*new_params.values())

                if new_error_rate < error_rate:
                    params = new_params
                    error_rate = new_error_rate
                    step_size[key] *= 1.2
                    continue

                # dan achteruit
                new_params[key] = params[key] - step_size[key]
                new_error_rate = self.mean_squared_error(*new_params.values())

                if new_error_rate < error_rate:
                    params = new_params
                    error_rate = new_error_rate
                    step_size[key] *= -1.2
                    continue

                # beide richtingen slechter → stap verkleinen
                step_size[key] *= 0.2

        return params


# ----------------------------------------------------------------------
# Demo / test
# ----------------------------------------------------------------------

def main():
    # "echte" parameters voor Von Bertalanffy
    real_states = [2.0, 1.0]

    # genereer "data" met Runge-Kutta
    solver = Solver(1000, 10, models.VonBertalanffyModel, 1e-7, 0, *real_states)
    time, real = solver.runge_kutta_function()

    # zoek parameters terug uit de data
    searcher = Searcher(real, time, 2, models.VonBertalanffyModel, predict_function="runge-kutta")
    found_vals = searcher.random_search_all()

    found_solver = Solver(1000, 10, models.VonBertalanffyModel, 1e-7, 0, *found_vals.values())
    f_time, found = found_solver.runge_kutta_function()

    plt.plot(time, real, "-k", label="Exact")
    plt.plot(f_time, found, ".:r", label="Random search fit")
    plt.axhline(0.0, lw=0.5, color="k")
    plt.axvline(0.0, lw=0.5, color="k")
    plt.grid(True)
    plt.legend()
    plt.xlabel("$t$")
    plt.ylabel("$V(t)$")
    plt.show()


if __name__ == "__main__":
    main()

