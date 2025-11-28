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

def surface_limited(V, c=1.0, d=10.0):
    return c * V / ((V + d) ** (1/3))

def allee(V, c=0.01, Vmin=10.0, Vmax=100.0):
    return c * (V - Vmin) * (Vmax - V)


def linear_limited(V, c=1.0, d=10.0):
    return c * V / (V + d)

def simulate(f, V0=1.0, t_end=50.0, n_steps=2000, **kwargs):
    t = 0.0
    V = V0
    dt = t_end / n_steps

    ts = [t]
    Vs = [V]

    for _ in range(n_steps):
        dV = f(V, **kwargs) * dt
        V = V + dV
        t = t + dt

        ts.append(t)
        Vs.append(V)

    return ts, Vs

def main():
    t_end = 100.0
    n_steps = 2000
    ts_allee, V_allee = simulate(
        allee,
        V0=20.0,
        t_end=t_end,
        n_steps=n_steps,
        c=0.001,
        Vmin=10.0,
        Vmax=100.0,
    )

    ts_lin, V_lin = simulate(
        linear_limited,
        V0=1.0,
        t_end=t_end,
        n_steps=n_steps,
        c=2.0,
        d=10.0,
    )

    ts_surf, V_surf = simulate(
        surface_limited,
        V0=1.0,
        t_end=t_end,
        n_steps=n_steps,
        c=1.0,
        d=10.0,
    )


    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    axes = axes.flatten()

    # Allee
    axes[0].plot(ts_allee, V_allee)
    axes[0].set_title("Allee-effect groei")
    axes[0].set_xlabel("t")
    axes[0].set_ylabel("V(t)")

    # Lineair gelimiteerd
    axes[1].plot(ts_lin, V_lin)
    axes[1].set_title("Lineair gelimiteerde groei")
    axes[1].set_xlabel("t")
    axes[1].set_ylabel("V(t)")

    # Oppervlakte-gelimiteerd
    axes[2].plot(ts_surf, V_surf)
    axes[2].set_title("Oppervlakte-gelimiteerde groei")
    axes[2].set_xlabel("t")
    axes[2].set_ylabel("V(t)")

    axes[3].plot(ts_allee, V_allee, label="Allee")
    axes[3].plot(ts_lin, V_lin, label="Lineair gelimiteerd")
    axes[3].plot(ts_surf, V_surf, label="Oppervlakte-gelimiteerd")
    axes[3].set_title("Alle modellen vergeleken")
    axes[3].set_xlabel("t")
    axes[3].set_ylabel("V(t)")
    axes[3].legend()

    fig.suptitle("Vergelijking van tumorgroeimodellen (alleen simulatie)", fontsize=14)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


if __name__ == "__main__":
    main()
