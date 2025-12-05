
from typing import Sequence, Dict, Optional
import numpy as np
from random import gauss
from .Solver import Solver
import numpy as np
import math

class Searcher:
    def __init__(
        self,
        real_vals: Sequence[float],
        time_values: Sequence[float],
        n_params: int,
        model,
        predict_function: Optional[str] = None,
    ):
        self.real_values = list(real_vals)
        self.time_values = list(time_values)
        self.n_params = n_params
        self.model = model
        self.predictor = predict_function

    def mean_squared_error(self, *params: float) -> float:
        try:
            solver = Solver(
                len(self.real_values) * 50,
                (self.time_values[-1] - self.time_values[0]),
                self.model,
                self.real_values[0],
                self.time_values[0],
                *params
            )

            if self.predictor == "runge-kutta":
                t_pred, V_pred = solver.runge_kutta_function()
            elif self.predictor == "heun":
                t_pred, V_pred = solver.heun_function()
            else:
                t_pred, V_pred = solver.euler_function()
        except Exception:
            return float("inf")

        # Interpolatie naar de echte tijdstippen
        try:
            V_interp = np.interp(self.time_values, t_pred, V_pred)
        except Exception:
            return float("inf")

        V_pred = np.asarray(V_interp, dtype=float)
        V_real = np.asarray(self.real_values, dtype=float)

        # sanity checks
        if np.any(~np.isfinite(V_pred)):
            return float("inf")

        # log-MSE, met clipping om log(0) te vermijden
        V_pred = np.clip(V_pred, 1e-9, None)
        V_real = np.clip(V_real, 1e-9, None)

        diff = np.log(V_pred) - np.log(V_real)
        mse = float(np.mean(diff * diff))
        return mse


    def random_search_all(
        self,
        sigma: float = 0.01,
        max_tries: int = 1000,
    ) -> Dict[int, float]:
        params = {i: 0.0 for i in range(self.n_params)}
        error_rate = self.mean_squared_error(*params.values())
        if not isinstance(error_rate, (int, float)):
            error_rate = float("inf")

        tries = 0
        while tries < max_tries:
            new_params = {
                key: val + gauss(0.0, sigma)
                for key, val in params.items()
            }
            new_error = self.mean_squared_error(*new_params.values())
            if not isinstance(new_error, (int, float)):
                new_error = float("inf")

            if new_error < error_rate:
                params = new_params
                error_rate = new_error
                tries = 0
            else:
                tries += 1

        return params