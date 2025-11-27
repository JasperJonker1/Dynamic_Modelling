from math import log
from matplotlib import pyplot as plt


def von_bertalanffy(V, t):
    d = 1
    c = 1
    return (c * pow(V, 2/3)) - (d * V)


def gompertz(V, c):
    cap = 1000
    c = 1
    return c * V * log(cap / V)


def main():
    t = 0
    volume = 1e-10

    ts = [t]
    ys = [volume]

    n_steps = 1000

    dt = 15 / n_steps

    for i in range(n_steps):
        dy1 = gompertz(volume, t) * dt
        t = t + dt
        volume = volume + dy1

        ts.append(t)
        ys.append(volume)

    plt.plot(ts, ys)
    plt.show()


if __name__ == "__main__":
    main()
