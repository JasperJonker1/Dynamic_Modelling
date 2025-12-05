from Solver import Solver
from Searcher import Searcher
import models
import sys
from math import log
import inspect
import csv


def main():
    """
    commandline input:
    python3 Main.py csv_file.csv search_function criterion_function [predict_function]
    search_function options: random / pattern
    criterion_function options: BIC / AIC / AICc
    predict_function [optional] options: runge-kutta / heun / euler (default)
    :return: prints mean squared error and the BIC, AIC or AICc (depending on input) of all models from given input
     to commandline.
    """
    user_args = sys.argv[1:]

    # get data from csv file
    with open(user_args[0]) as csv_file:
        volume, time_steps = csv.reader(csv_file)

    # convert to the right data types with list comprehension
    volume = [float(i) for i in volume]
    time_steps = [int(i) for i in time_steps]

    n_steps = len(volume)
    time = time_steps[-1] - time_steps[0]

    # define all models
    all_models = [models.AlleeModel, models.SurfaceLimitedModel, models.VonBertalanffyModel, models.GompertzPaperModel,
                  models.GompertzLesModel, models.LinearLimitedModel]

    mean_squared_error = {}
    information_criterion = {}

    predict_function = None
    # assign predict function if given
    if len(user_args) > 3:
        predict_function = user_args[3]

    # loop through all models and add the MSE and information criterion to the dicts
    for model in all_models:
        k = len(inspect.signature(model.__init__).parameters) - 1
        searcher = Searcher(volume, time_steps, k, model, predict_function)
        if user_args[1] == 'random':
            params = searcher.random_search_all()
        elif user_args[1] == 'pattern':
            params = searcher.pattern_search_all()
        else:
            raise SyntaxError(f"Argument: {user_args[1]} is invallid")

        predicted_solver = Solver(n_steps, time, model, volume[0], time_steps[0], *params)
        predicted_time, predicted_volume = predicted_solver.runge_kutta_function()

        mean_squared_error[model.__name__] = searcher.mean_squared_error(*params)

        if user_args[2] == "BIC":
            information_criterion[model.__name__] = (n_steps * log(searcher.mean_squared_error(*params)) +
                                                     k * log(len(predicted_volume)))
        elif user_args[2] == "AIC":
            information_criterion[model.__name__] = (n_steps * log(searcher.mean_squared_error(*params)) +
                                                     k * 2)
        elif user_args[2] == "AICc":
            information_criterion[model.__name__] = (n_steps * log(searcher.mean_squared_error(*params)) +
                                                     k * 2 * (n_steps / (n_steps - k - 1)))
    print(mean_squared_error)
    print(information_criterion)


if __name__ == "__main__":
    main()
