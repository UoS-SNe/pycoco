"""

"""
import numpy as np


__all__ = ["bazin09", "karpenka12", "firth17",
           "bazin09_listarg", "karpenka12_listarg", "firth17_listarg",
           "_defined_models"]


_defined_models = ["bazin09", "karpenka12", "firth17"]


def bazin09(x, a, t_0, t_rise, t_fall):
    return a * np.exp(-(x - t_0) / t_fall) / (1.0 + np.exp(- (x - t_0) / t_rise))


def bazin09_listarg(x, params):
    return params[0] * np.exp(-(x - params[1]) / params[3]) / (1.0 + np.exp(- (x - params[1]) / params[2]))


def karpenka12(x, a, b, t_0, t_1, t_rise, t_fall):
    return a * (1. + b * (x - t_1)*(x - t_1)) * np.exp(- (x - t_0)/t_fall) / (1. + np.exp(- (x - t_0) / t_rise))


def karpenka12_listarg(x, params):
    return params[0] * (1. + params[1] * (x - params[3])*(x - params[3])) * np.exp(- (x - params[2])/params[5]) / (1. + np.exp(- (x - params[2]) / params[4]))


def firth17(x, a, b, t_0, t_1, t_2, t_x, t_rise, t_fall):
    numerator = a * (1. + b * (x - t_1) * (x - t_1)) * np.exp(- (x - t_0) / t_fall)
    denominator = (1. + np.exp(- (x - t_0) / t_rise)) / (1. + np.exp(- (x - t_2) / t_x))
    return numerator/denominator


def firth17_listarg(x, params):
    numerator = params[0] * (1. + params[1] * (x - params[3]) * (x - params[3])) * np.exp(- (x - params[2]) / params[7])
    denominator = (1. + np.exp(- (x - params[2]) / params[6])) / (1. + np.exp(- (x - params[4]) / params[5]))
    return numerator/denominator

