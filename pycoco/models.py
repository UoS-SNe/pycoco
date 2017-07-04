import numpy as np


__all__ = ["bazin09", "karpenka12", "firth17"]


def bazin09(x, a, t_0, t_rise, t_fall):
    return a * np.exp(-(x - t_0) / t_fall) / (1.0 + np.exp(- (x - t_0) / t_rise))


def karpenka12(x, a, b, t_0, t_1, t_rise, t_fall):
    return a * (1. + b * (x - t_1)*(x - t_1)) * np.exp(- (x - t_0)/t_fall) / (1. + np.exp(- (x - t_0) / t_rise))


def firth17(x, a, b, t_0, t_1, t_2, t_x, t_rise, t_fall):
    numerator = a * (1. + b * (x - t_1) * (x - t_1)) * np.exp(- (x - t_0) / t_fall)
    denominator = (1. + np.exp(- (x - t_0) / t_rise)) / (1. + np.exp(- (x - t_2) / t_x))
    return numerator/denominator
