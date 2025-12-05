import math
from math import log
from abc import ABC, abstractmethod


class TumorGrowthModel(ABC):
    """Abstract basis voor alle tumorgroeimodellen."""

    @abstractmethod
    def dVdt(self, V: float, t: float) -> float:
        pass


class LinearGrowthModel(TumorGrowthModel):
    """dV/dt = c"""
    def __init__(self, c=1.0):
        self.c = c

    def dVdt(self, V, t):
        return self.c


class ExponentialGrowthModel(TumorGrowthModel):
    """dV/dt = c * V"""
    def __init__(self, c=1.0):
        self.c = c

    def dVdt(self, V, t):
        return self.c * V


class MendelsohnGrowthModel(TumorGrowthModel):
    """dV/dt = c * V^d"""
    def __init__(self, c=1.0, d=1.0):
        self.c = c
        self.d = d

    def dVdt(self, V, t):
        if V < 0:
            V = 0.0
        if V > 1e6:
            V = 1e6
        return self.c * (V ** self.d)


class ExponentialSaturatingModel(TumorGrowthModel):
    """dV/dt = c * (Vmax - V)"""
    def __init__(self, c=1.0, Vmax=1.0):
        self.c = c
        self.Vmax = Vmax

    def dVdt(self, V, t):
        return self.c * (self.Vmax - V)


class LogisticGrowthModel(TumorGrowthModel):
    """dV/dt = c * V * (Vmax - V)"""
    def __init__(self, c=1.0, Vmax=1.0):
        self.c = c
        self.Vmax = Vmax

    def dVdt(self, V, t):
        return self.c * V * (self.Vmax - V)


class MontrollGrowthModel(TumorGrowthModel):
    """dV/dt = c * V * (Vmax^d - V^d)"""
    def __init__(self, c=1.0, Vmax=1.0, d=1.0):
        self.c = c
        self.Vmax = Vmax
        self.d = d

    def dVdt(self, V, t):
        return self.c * V * (self.Vmax ** self.d - V ** self.d)


class VonBertalanffyModel(TumorGrowthModel):
    """dV/dt = c * V^(2/3) - d * V"""
    def __init__(self, c=1.0, d=1.0):
        self.c = c
        self.d = d

    def dVdt(self, V, t):
        if V < 0:
            V = 0.0
        return self.c * (V ** (2/3)) - self.d * V


class GompertzLesModel(TumorGrowthModel):
    """dV/dt = c * V * ln(cap / V)"""
    def __init__(self, c=1.0, cap=1000):
        self.c = c
        self.cap = cap if cap > 0 else 1e-6

    def dVdt(self, V, t):
        if V <= 0:
            V = 1e-9
        return self.c * V * log(self.cap / V)


class GompertzPaperModel(TumorGrowthModel):
    """dV/dt = alpha * exp(-beta t) * V"""
    def __init__(self, alpha=1.0, beta=1.0):
        self.alpha = alpha
        self.beta = beta

    def dVdt(self, V, t):
        if V < 0:
            V = 0.0
        return self.alpha * math.exp(-self.beta * t) * V


class SurfaceLimitedModel(TumorGrowthModel):
    """dV/dt = c * V / (V + d)^(1/3)"""
    def __init__(self, c=1.0, d=10.0):
        self.c = c
        self.d = d

    def dVdt(self, V, t):
        return self.c * V / ((V + self.d) ** (1/3))


class AlleeModel(TumorGrowthModel):
    """dV/dt = c * (V - Vmin) * (Vmax - V)"""
    def __init__(self, c=0.01, Vmin=10.0, Vmax=100.0):
        self.c = c
        self.Vmin = Vmin
        self.Vmax = Vmax

    def dVdt(self, V, t):
        return self.c * (V - self.Vmin) * (self.Vmax - V)


class LinearLimitedModel(TumorGrowthModel):
    """dV/dt = c * V / (V + d)"""
    def __init__(self, c=1.0, d=10.0):
        self.c = c
        self.d = d

    def dVdt(self, V, t):
        return self.c * V / (V + self.d)
