from tumor_models import Solver, Searcher
from tumor_models import (
    AlleeModel,
    SurfaceLimitedModel,
    VonBertalanffyModel,
    GompertzPaperModel,
    GompertzLesModel,
    LinearLimitedModel,
)
import sys
import csv
import inspect
from math import log


def main():
    """
    Gebruik:
        python3 Main.py data.csv search_function criterion [predict_function]

    Waar:
      - data.csv : CSV met twee regels:
           rij 1: volumes  (V1,V2,...) 
           rij 2: tijdstippen (t1,t2,...)
      - search_function : random  of  pattern
      - criterion       : BIC / AIC / AICc
      - predict_function (optioneel): runge-kutta / heun / euler (default: euler)

    Output:
      - dictionary met MSE per model
      - dictionary met gekozen informatiecriterium per model
    """

    user_args = sys.argv[1:]
    if len(user_args) < 3:
        print("Gebruik: Main.py data.csv search_function (random|pattern) criterion (BIC|AIC|AICc) [predict_function]")
        sys.exit(1)

    csv_path = user_args[0]
    search_method = user_args[1]
    criterion = user_args[2]      
    predict_function = user_args[3] if len(user_args) > 3 else None


    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)

    if len(rows) < 2:
        raise ValueError("CSV moet minstens twee rijen bevatten: eerste volumes, tweede tijdstippen.")

    # hier ga ik uit van: eerste rij = volumes, tweede rij = tijdstippen
    volume = [float(x) for x in rows[0] if x.strip() != ""]
    time_steps = [float(x) for x in rows[1] if x.strip() != ""]

    if len(volume) != len(time_steps):
        raise ValueError("Aantal volumes en tijdstippen in CSV komt niet overeen.")

    n = len(volume)
    total_time = time_steps[-1] - time_steps[0]


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

        searcher = Searcher(
            real_vals=volume,
            time_values=time_steps,
            n_params=k,
            model=model,
            predict_function=predict_function,  # None/euler, heun, runge-kutta
        )

        # parameters zoeken
        if search_method == "random":
            param_dict = searcher.random_search_all()
        elif search_method == "pattern":
            param_dict = searcher.pattern_search_all()
        else:
            raise SyntaxError(f"Argument: {search_method} is ongeldig (gebruik 'random' of 'pattern').")

        params = list(param_dict.values())

        # MSE op de data
        mse = searcher.mean_squared_error(*params)
        mse_dict[model.__name__] = mse

        # informatiecriterium: gebruik n = aantal observaties
        if criterion == "BIC":
            ic = n * log(mse) + k * log(n)
        elif criterion == "AIC":
            ic = n * log(mse) + 2 * k
        elif criterion == "AICc":
            ic = n * log(mse) + 2 * k * (n / (n - k - 1))
        else:
            raise ValueError(f"Onbekend criterium: {criterion}")

        ic_dict[model.__name__] = ic

    print("MSE per model:")
    for name, mse in mse_dict.items():
        print(f"  {name:20s}: {mse:.6g}")

    print(f"\n{criterion} per model:")
    for name, ic in ic_dict.items():
        print(f"  {name:20s}: {ic:.6g}")


if __name__ == "__main__":
    main()
