"""

"""

import numpy as np

__all__ = ["bazin09",
           "karpenka12",
           "firth17",
           ]

def bazin09(X, A, t_0, t_rise, t_fall):
    """

    """

    return A * np.exp(-(X-t_0)/t_fall) / (1.0 + np.exp(- (X - t_0)/t_rise))

def karpenka12(X, A, B, t_0, t_1, t_rise, t_fall):
    """

    """

    return A * (1. + B * (X - t_1)*(X - t_1)) * np.exp(- (X - t_0)/t_fall) / (1. + np.exp(- (X - t_0) / t_rise ))

def firth17(X, A, B, t_0, t_1, t_2, t_x, t_rise, t_fall):
    """

    """

    return A * (1. + B * (X - t_1)*(X - t_1)) * np.exp(- (X - t_0)/t_fall) / (1. + np.exp(- (X - t_0) / t_rise )) / (1. + np.exp(- (X - t_2) / t_x ))
