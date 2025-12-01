from random import gauss
import Algorithm_calculations


class Searcher:

    def __init__(self, real_vals, n_params, predict_function, algorithm):
        self.real_values = real_vals
        self.n_params = n_params
        self.predictor = predict_function
        self.algorithm = algorithm

    def mean_squared_error(self, **params):
        ts, predicted_values = self.predictor(params, self.algorithm)
        mse = sum((predict_val - real_val) ** 2 for predict_val, real_val in zip(predicted_values, self.real_values))
        mse = mse / (len(ts))
        return mse

    def random_search(self):
        """
        finds optimal parameters of a curve by randomly changing the parameters
        :return: found optimal parameters
        """
        params = {i: 0.0 for i in range(self.n_params)}
        error_rate = self.mean_squared_error(**params)

        # loop through the parameters to find the best of each individual parameter
        for i in range(self.n_params):
            tries = 0
            while True:
                tries += 1
                new_params = params.copy()
                # create new parameter randomly
                new_params[i] = params[i] + gauss(mu=1.0, sigma=1.0)
                new_error_rate = self.mean_squared_error(**new_params)
                # replace parameter if new parameter gives lower error rate
                if new_error_rate < error_rate:
                    params[i] = new_params[i]
                    error_rate = new_error_rate
                    tries = 0
                if error_rate < 1e-6 or tries > 10000:
                    break
        return params

    def pattern_search(self):
        """
        finds optimal parameters by changing parameters
        :return: calculated optimal parameters
        """
        params = {i: 0.0 for i in range(self.n_params)}
        error_rate = self.mean_squared_error(**params)

        for i in range(self.n_params):
            step_size = 10.0
            while True:
                new_params = params.copy()
                new_params[i] = params[i] + step_size
                new_error_rate = self.mean_squared_error(**new_params)
                if new_error_rate < error_rate:
                    params[i] = new_params[i]
                    error_rate = error_rate
                else:
                    step_size *= -0.9
                if error_rate < 1e-6 or -1e-2 < step_size > 1e-2:
                    break
        return params


def main():
    print("Hello World!")


if __name__ == "__main__":
    main()
