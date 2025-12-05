"""
tumor_models package

Bevat:
- groeimodellen (Gompertz, Von Bertalanffy, Logistic, enz.)
- ODE-oplossers (Euler, Heun, Runge-Kutta) via Solver
- Fit-functionaliteit (MSE, random search) via Searcher
"""

# === Models ===
from .models import (
    TumorGrowthModel,
    VonBertalanffyModel,
    GompertzLesModel,
    GompertzPaperModel,
    SurfaceLimitedModel,
    AlleeModel,
    LinearLimitedModel,
    LinearGrowthModel,
    ExponentialGrowthModel,
    MendelsohnGrowthModel,
    ExponentialSaturatingModel,
    LogisticGrowthModel,
    MontrollGrowthModel,
)

# === Solver ===
from .Solver import Solver

# === Searcher (fitting) ===
from .Searcher import Searcher


__all__ = [
    # Models
    "TumorGrowthModel",
    "VonBertalanffyModel",
    "GompertzLesModel",
    "GompertzPaperModel",
    "SurfaceLimitedModel",
    "AlleeModel",
    "LinearLimitedModel",
    "LinearGrowthModel",
    "ExponentialGrowthModel",
    "MendelsohnGrowthModel",
    "ExponentialSaturatingModel",
    "LogisticGrowthModel",
    "MontrollGrowthModel",

    # Solver + Searcher
    "Solver",
    "Searcher",
]
