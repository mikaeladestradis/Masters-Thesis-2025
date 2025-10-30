from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Any, List, Callable
import itertools
import math

# find the fold points, defined as the points at which the partial 
# derivative of the input function with respect to the given variable is 0.
# this function works best in the two dimensional case
def find_fold_points(function, graph, w, variable, params, guesses, stimulus):
    fold_points = []
    func = diff(function, variable)
    func = func.subs(w, graph)
    func = func.subs(params)
    func = func.subs('I', stimulus)
    for guess in guesses: 
        fold_point = float(nsolve(func, variable, guess))
        if fold_point not in fold_points:
            fold_points.append(fold_point)
    return fold_points

def find_fold_lines(fold_function, v_vals, z_vals, variables, stimulus, params):
    v = variables[0]
    z = variables[1]

    # given the function, determine all the points at which the value is zero.
    fold = fold_function.subs('I', stimulus)
    fold = fold.subs(params)
    fold_surface = lambdify((v,z), fold)
    
    v_fold = []
    z_fold = []
    for i, j in itertools.product(v_vals.tolist(), z_vals.tolist()):
        if math.isclose(fold_surface(i,j), 0, abs_tol = 0.0011):
            v_fold.append(i)
            z_fold.append(j)
            if j == z_vals[-1]:
                z_fold[-1] = np.nan
                v_fold[-1] = np.nan
    
    return [v_fold, z_fold]

#as input: a function of only ONE variable, and a list of parameter values for all the other variables
def find_equilibria(func, var, guesses: float, stimulus: float) -> List:
    equilibria = []
    func = func.subs('I', stimulus)
    for guess in guesses:
        try: 
            next_eq = float(nsolve(func, var, guess))
            if next_eq not in equilibria:
                equilibria.append(next_eq)
        except:
            pass
    return equilibria

#lambdify versions;
def convert_funcs(functions: List[Any], vars):
    converted = []
    for func in functions:
        converted.append(lambdify(vars, func, modules='numpy'))
    return converted

#Plotting!
def plot(functions: List[tuple], add_folds:bool, folds:List[float], label):
    x_vals = np.linspace(-0.99, 0.6, 1000)
    plt.figure(figsize=(10, 6))

    for (func, params) in functions:
        try:
            w_vals = func(x_vals, **params)
        except:
            w_vals = func(x_vals)
        plt.plot(w_vals, x_vals, color='black')

    if add_folds:
        plt.plot(func(folds[0], **params), folds[0], marker='o', color='red', label='F1')
        plt.plot(func(folds[1], **params), folds[1], marker='o', color='blue', label='F2')
    plt.xlabel('w')
    plt.ylabel('v')
    plt.title(label)
    plt.grid()
    plt.xlim(-0.2, 0.4)
    plt.legend()
    plt.show()
    plt.savefig(label)

#translation of functions, this we need to use the normal functions!
def transformation(functions: List[Any], transformation: Any, variable: Any) -> List[Any]:
    transformed_functions = []
    for func in functions: 
        transformed_functions.append(func.subs(variable, transformation))
    return transformed_functions
