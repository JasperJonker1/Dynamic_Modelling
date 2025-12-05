# cli_compare_models.py
from tumor_models import Solver, Searcher
from tumor_models import (
    AlleeModel,
    SurfaceLimitedModel,
    VonBertalanffyModel,
    GompertzPaperModel,
    GompertzLesModel,
    LinearLimitedModel,
)
import math
import inspect
import sys


def main():
    """
    Gebruik:
        python cli_compare_models.py <n_steps> <time> <init_volume> random <BIC|AIC|AICc>

    Voorbeeld:
        python cli_compare_models.py 1000 10 1e-7 random BIC
    """
    user_args = sys.argv[1:]
    if len(user_args) != 5:
        print("Gebruik: n_steps time init_volume random (BIC|AIC|AICc)")
        sys.exit(1)

    n_steps = int(user_args[0])
    total_time = float(user_args[1])
    init_volume = float(user_args[2])
    method = user_args[3]
    crit = user_args[4]

    # synthetic "data" met Von Bertalanffy
    init_time = 0.0
    true_params = (2.0, 1.6)
    data_solver = Solver(n_steps, total_time, VonBertalanffyModel,
                         init_volume, init_time, *true_params)
    time_steps, volume = data_solver.runge_kutta_function()

    all_models = [
        AlleeModel,
        SurfaceLimitedModel,
        VonBertalanffyModel,
        GompertzPaperModel,
        GompertzLesModel,
        LinearLimitedModel,
    ]

    mse_dict = {}
    ic_dict = {}

    for model in all_models:
        k = len(inspect.signature(model.__init__).parameters) - 1

        searcher = Searcher(volume, time_steps, k, model, predict_function="runge-kutta")

        if method == "random":
            params = searcher.random_search_all()
            param_list = list(params.values())
        else:
            raise SyntaxError(f"Argument: {method} is ongeldig")

        mse = searcher.mean_squared_error(*param_list)
        mse_dict[model.__name__] = mse

        n = len(volume)

        if crit == "BIC":
            ic = n * math.log(mse) + k * math.log(n)
        elif crit == "AIC":
            ic = n * math.log(mse) + 2 * k
        elif crit == "AICc":
            ic = n * math.log(mse) + 2 * k * (n / (n - k - 1))
        else:
            raise ValueError(f"Onbekend criterium: {crit}")

        ic_dict[model.__name__] = ic

    print("MSE per model:")
    for name, mse in mse_dict.items():
        print(f"  {name:20s}  {mse:.6g}")

    print(f"\n{crit} per model:")
    for name, ic in ic_dict.items():
        print(f"  {name:20s}  {ic:.6g}")


if __name__ == "__main__":
    main()

