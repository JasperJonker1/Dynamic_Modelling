from Solver import Solver
from Searcher import Searcher
import models
import sys
from math import log
import inspect


def main():
    """
    :return: prints mean squared error and the BIC, AIC or AICc (depending on input) of all models from given input
     to commandline.
    """
    user_args = sys.argv[1:]
    print(user_args)
    n_steps = int(user_args[0])
    time = int(user_args[1])
    init_volume = float(user_args[2])
    init_time = 0
    solver = Solver(n_steps, time, models.VonBertalanffyModel, init_volume, init_time, 2.0, 1.6)
    time_steps, volume = solver.runge_kutta_function()

    all_models = [models.AlleeModel, models.SurfaceLimitedModel, models.VonBertalanffyModel, models.GompertzPaperModel,
                  models.GompertzLesModel, models.LinearLimitedModel]

    mean_squared_error = {}
    information_criterion = {}

    for model in all_models:
        k = len(inspect.signature(model.__init__).parameters) - 1
        searcher = Searcher(volume, time_steps, k, model)
        if user_args[3] == 'random':
            params = searcher.random_search_all()
        else:
            raise SyntaxError(f"Argument: {user_args[3]} is invallid")

        predicted_solver = Solver(n_steps, time, model, init_volume, init_time, *params)
        predicted_time, predicted_volume = predicted_solver.runge_kutta_function()

        mean_squared_error[model.__name__] = searcher.mean_squared_error(*params)

        if user_args[4] == "BIC":
            information_criterion[model.__name__] = (n_steps * log(searcher.mean_squared_error(*params)) +
                                                     k * log(len(predicted_volume)))
        elif user_args[4] == "AIC":
            information_criterion[model.__name__] = (n_steps * log(searcher.mean_squared_error(*params)) +
                                                     k * 2)
        elif user_args[4] == "AICc":
            information_criterion[model.__name__] = (n_steps * log(searcher.mean_squared_error(*params)) +
                                                     k * 2 * (n_steps / (n_steps - k - 1)))
    print(mean_squared_error)
    print(information_criterion)


if __name__ == "__main__":
    main()
