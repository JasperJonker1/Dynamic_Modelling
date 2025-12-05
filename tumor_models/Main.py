from Solver import Solver
from Searcher import Searcher
import models
import sys
from math import log, isfinite
import inspect
import csv


def main():
    """
    Gebruik:
      python3 Main.py data.csv search_function criterion [predict_function]

    search_function : random | pattern
    criterion       : BIC | AIC | AICc
    predict_function (optioneel) : runge-kutta | heun | euler (default: euler/None)
    """

    user_args = sys.argv[1:]
    if len(user_args) < 3:
        print("Gebruik: Main.py data.csv search_function (random|pattern) criterion (BIC|AIC|AICc) [predict_function]")
        sys.exit(1)

    csv_path = user_args[0]
    search_method = user_args[1]
    criterion = user_args[2]
    predict_function = user_args[3] if len(user_args) > 3 else None

    # --- CSV inlezen ---
    with open(csv_path, newline="") as csv_file:
        reader = csv.reader(csv_file)
        rows = list(reader)

    if len(rows) < 2:
        raise ValueError("CSV moet minstens twee rijen bevatten: eerste rij volumes, tweede rij tijdstippen.")

    # eerste rij = volumes, tweede rij = tijd
    volume = [float(i) for i in rows[0] if i.strip() != ""]
    time_steps = [float(i) for i in rows[1] if i.strip() != ""]

    if len(volume) != len(time_steps):
        raise ValueError("Aantal volumes en tijdstippen in CSV komt niet overeen.")

    n_obs = len(volume)
    total_time = time_steps[-1] - time_steps[0]

    # --- modellenlijst ---
    all_models = [
        models.AlleeModel,
        models.SurfaceLimitedModel,
        models.VonBertalanffyModel,
        models.GompertzPaperModel,
        models.GompertzLesModel,
        models.LinearLimitedModel,
    ]

    mse_dict = {}
    ic_dict = {}

    # --- modellen fitten ---
    for model in all_models:
        # aantal parameters = #__init__-argumenten - 1 (self)
        k = len(inspect.signature(model.__init__).parameters) - 1

        searcher = Searcher(
            real_vals=volume,
            time_values=time_steps,
            n_params=k,
            model=model,
            predict_function=predict_function,
        )

        # parameters zoeken
        if search_method == "random":
            param_dict = searcher.random_search_all()
        elif search_method == "pattern":
            # alleen gebruiken als je pattern_search_all geïmplementeerd hebt
            param_dict = searcher.pattern_search_all()
        else:
            raise SyntaxError(f"Argument: {search_method} is ongeldig (gebruik 'random' of 'pattern').")

        params = list(param_dict.values())

        # MSE op de data
        mse = searcher.mean_squared_error(*params)
        mse_dict[model.__name__] = mse

        # informatiecriterium (alleen als MSE zinvol is)
        if (not isfinite(mse)) or mse <= 0:
            ic_dict[model.__name__] = float("inf")
            continue

        if criterion == "BIC":
            ic = n_obs * log(mse) + k * log(n_obs)
        elif criterion == "AIC":
            ic = n_obs * log(mse) + 2 * k
        elif criterion == "AICc":
            ic = n_obs * log(mse) + 2 * k * (n_obs / (n_obs - k - 1))
        else:
            raise ValueError(f"Onbekend criterium: {criterion}")

        ic_dict[model.__name__] = ic

    # --- resultaten printen ---
    print("MSE per model:")
    for name, mse in mse_dict.items():
        print(f"  {name:20s}: {mse:.6g}")

    print(f"\n{criterion} per model:")
    for name, ic in ic_dict.items():
        print(f"  {name:20s}: {ic:.6g}")


if __name__ == "__main__":
    main()
