from random import gauss
from Solver import Solver
import models
from matplotlib import pyplot as plt


class Searcher:

    def __init__(self, real_vals, time_values, n_params, model, predict_function=None):
        self.real_values = real_vals
        self.time_values = time_values
        self.n_params = n_params
        self.model = model
        self.predictor = predict_function

    def mean_squared_error(self, *params):
        solver = Solver(len(self.real_values) + 1, (self.time_values[-1] - self.time_values[0]), self.model,
                        self.real_values[0], self.time_values[0], *params)
        if self.predictor == "runge-kutta":
            ts, predicted_values = solver.runge_kutta_function()
        elif self.predictor == "heun":
            ts, predicted_values = solver.heun_function()
        else:
            ts, predicted_values = solver.euler_function()

        sum_squared_error = 0.0
        for predict_val, real_val in zip(predicted_values, self.real_values):
            error = predict_val - real_val
            sum_squared_error += error * error
        mse = sum_squared_error / (len(predicted_values))
        return abs(mse)

    def random_search_one(self):
        """
        finds optimal parameters of a curve by randomly changing the parameters
        :return: found optimal parameters
        """
        params = {i: 0.0 for i in range(self.n_params)}
        error_rate = self.mean_squared_error(*params)

        # loop through the parameters to find the best of each individual parameter
        for i in range(self.n_params):
            tries = 0
            while True:
                tries += 1
                new_params = params.copy()
                # create new parameter randomly
                new_params[i] = params[i] + gauss(mu=1.0, sigma=1.0)
                new_error_rate = self.mean_squared_error(*new_params.values())
                # replace parameter if new parameter gives lower error rate
                if new_error_rate < error_rate:
                    params[i] = new_params[i]
                    error_rate = new_error_rate
                    tries = 0
                if tries > 10000:
                    break
        return params

    def pattern_search(self):
        """
        finds optimal parameters by changing parameters
        :return: calculated optimal parameters
        """
        params = {i: 0.0 for i in range(self.n_params)}
        error_rate = self.mean_squared_error(*params)

        for i in range(self.n_params):
            step_size = 10.0
            while True:
                new_params = params.copy()
                new_params[i] = params[i] + step_size
                new_error_rate = self.mean_squared_error(*new_params.values())
                if new_error_rate < error_rate:
                    params[i] = new_params[i]
                    error_rate = new_error_rate
                else:
                    step_size *= -0.8
                if -1e-6 < step_size < 1e-6:
                    break
        return params

    def random_search_all(self):
        """
        finds optimal parameters of a curve by randomly changing the parameters
        :return: found optimal parameters
        """
        params = {i: 0.0 for i in range(self.n_params)}
        error_rate = self.mean_squared_error(*params)

        tries = 0
        while tries < 1000:
            new_params = {key: val + gauss(sigma=0.01) for key, val in params.items()}
            new_error_rate = self.mean_squared_error(*new_params.values())
            if new_error_rate < error_rate:
                params = new_params
                error_rate = new_error_rate
                tries = 0
            tries += 1
        return params

    def pattern_search_all(self):
        """
        finds optimal parameters by changing parameters
        :return: calculated optimal parameters
        """
        # set default
        params = {i: 0.0 for i in range(self.n_params)}
        step_size = {i: 10.0 for i in range(self.n_params)}
        error_rate = self.mean_squared_error(*params)

        while max(step_size.values()) < 1e-6 or min(step_size.values()) > -1e-6:
            for key in params:
                new_params = params.copy()
                new_params[key] = params[key] + step_size[key]
                new_error_rate = self.mean_squared_error(*new_params.values())
                if new_error_rate < error_rate:
                    params = new_params
                    error_rate = new_error_rate
                    step_size[key] *= 1.2
                    continue
                new_params[key] = params[key] - step_size[key]
                new_error_rate = self.mean_squared_error(*new_params.values())
                if new_error_rate < error_rate:
                    params = new_params
                    error_rate = new_error_rate
                    step_size[key] *= -1.2
                    continue
                step_size[key] *= 0.2
        return params


def main():
    real_states = [2.0, 1.0]
    solver = Solver(1000, 10, models.VonBertalanffyModel, 1e-7, 0, *real_states)
    time, real = solver.runge_kutta_function()
    # print(real)

    searcher = Searcher(real, time, 2, models.VonBertalanffyModel)
    found_vals = searcher.random_search_all()

    found_solver = Solver(1000, 10, models.VonBertalanffyModel, 1e-7, 0, *found_vals.values())
    f_time, found = found_solver.runge_kutta_function()
    # print(found)

    # print(found_vals)
    plt.plot(time, real, '-k', label='Exact')
    plt.plot(time[-1], real[-1], 'ok')
    plt.plot(f_time, found, '.:r', label='Direct search')
    plt.axhline(0.0, lw=0.5, color='k')
    plt.axvline(0.0, lw=0.5, color='k')
    plt.grid(True)
    plt.legend()
    plt.xlabel('$t$')
    plt.ylabel('$y(t)$')
    plt.show()


if __name__ == "__main__":
    main()
