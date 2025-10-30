import sympy as sp
import taylor as T
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d.axes3d import Axes3D
import util
import warnings
from scipy.integrate import solve_ivp
from scipy.differentiate import jacobian
from systems.system import System
import matplotlib.style as mplstyle
# ok so actually, I need to fix this as it is wrong!!!

# tonight: I want to get the folds plotted properly, as in maybe I can take a function of z but for split v values 
# so two functions almost idk if I can do that...

warnings.filterwarnings('ignore')
points = 1000
kv = 100
int_method = "LSODA"

#for presentation - no transient (~20), one transient (~50), two transient (~80) three? #82.355
I_stim = 82.35/(kv*20) 
eq_guesses = [-80/kv, -40/kv, -50/kv, -10/kv]
fold_guesses = [10/kv, -60/kv]

# I really really need to go through and check something, it seems I can only get -y when I have positive z,and 
# negative z when I have positive y, which should NOT be the case?
# V = 100v (where 100 is mV)
y, z, v = sp.symbols('y z v')
I, N, a, m, K, a_k, sub, S, a_l, L = sp.symbols('I, N, a, m, K, a_k, sub, S, a_l, L')
eps, gamma_y, beta_y, gamma_z, beta_z = sp.symbols('epsilon, gamma_y, beta_y, gamma_z, beta_z')
dvdt, dydt, dzdt, yinf, zinf = sp.symbols('dv/dt, dy/dt, dz/dt, yinf, zinf')

yinf = 1/2*(1+sp.tanh(gamma_y*v - beta_y))
zinf = 1/2*(1+sp.tanh(gamma_z*v - beta_z))
dvdt = I - 1/2*(1+sp.tanh(a*v-m))*(v - N) - a_k*y*(v-K) - sub*z*(v-S) - a_l*(v-L)
dydt = eps*(yinf - y)*sp.cosh((gamma_y*v - beta_y)/2)
dzdt = eps*(zinf - z)*sp.cosh((gamma_z*v - beta_z)/2)
phi = (I - 1/2*(1+sp.tanh(a*v-m))*(v - N) - sub*z*(v-S) - a_l*(v-L))/(a_k*(v-K))

#todo:  need to edit this so it is in the non-dimensionalised variables!!!!
params = {
    'I': I_stim,
    'epsilon': 0.00015, # then this will change to C/gf * phi_z = 0.15* (2/20) ~ 0.015
    'N': 50/kv,
    'a': kv/18,
    'm': -1.2/18,
    'K': -100/kv,
    'a_k': 20/20,
    'sub': 8/20,
    'S': -100/kv,
    'a_l': 2.0/20,
    'L': -70/kv,
    'gamma_y': kv/10,
    'beta_y': -10/10, #now we are changing this one...
    'gamma_z': kv/15,
    'beta_z': -21/15,
}

# different parameter setups. 
sympy_params = {sp.symbols(k): v for k, v in params.items()}
sympy_params_no_I = {sp.symbols(k): v for k, v in params.items() if k != 'I'}
sympy_params_no_eps = {sp.symbols(k): v for k,v in params.items() if k != 'epsilon'}

# this is to compute and plot the function y = phi
phi_subs = phi.subs(sympy_params)
phi_subs_no_I = phi.subs(sympy_params_no_I)
plot_phi = sp.lambdify((v, z), phi_subs)

# to plot the y 'nullcline'
yinf_subs = yinf.subs(sympy_params)
yinf_subs_og = yinf.subs(sympy_params_no_I)
yinf_subs_og = yinf_subs_og.subs('I', 0.0)
plot_yinf = sp.lambdify(v, yinf_subs)
plot_yinf_og = sp.lambdify(v, yinf_subs_og)

# to plot the z 'nullcline'
zinf_subs = zinf.subs(sympy_params)
zinf_subs_og = zinf.subs(sympy_params_no_I)
zinf_subs_og = zinf_subs_og.subs("I", 0.0)
plot_zinf = sp.lambdify(v, zinf_subs)
plot_zinf_og = sp.lambdify(v, zinf_subs_og)

# equilibria occur when: f = g = h = 0: ie) when  y =YINF(v), z=ZINF(v) and f=0
# original or 'og' equilibria are compute when the stimulus is 0, and this is our starting point. 
eq_curve = dvdt.subs(y, yinf)
eq_curve = eq_curve.subs(z, zinf)
eq_curve_subs = eq_curve.subs(sympy_params_no_I)
og_equilibria = util.find_equilibria(eq_curve_subs, v, eq_guesses, 0.0)
equilibria = util.find_equilibria(eq_curve_subs, v, eq_guesses, I_stim)

# finding the fold lines: occur when (dvdt=f) -> df/dv = 0.
# the fold lines are stored as a list, with the first index a list of v values 
# and the second a list of z values.
z_vals = np.linspace(-0.3, 0.3, points) 
v_vals = np.linspace(-0.9, 0.4, points)

d_v_f = sp.diff(dvdt, v)
eq_fold = d_v_f.subs(y, phi)
fold_lines = util.find_fold_lines(eq_fold, v_vals, z_vals, [v, z], I_stim, sympy_params_no_I)

# Also want the folded singularity, which occurs when the Dvf = 0, and
# Dyf * h + Dzf * g = 0 curves intersect. 
FSN_values = []

d_y_f = sp.diff(dvdt, y)
d_z_f = sp.diff(dvdt, z)
FSN_eq = d_y_f*(dydt) + d_z_f*(dzdt)
FSN_eq = FSN_eq.subs(y, phi)
FSN_eq_subs = FSN_eq.subs(sympy_params_no_I)
FSN_eq_with_stim = FSN_eq_subs.subs(I, I_stim)
reduced_nullcline = sp.lambdify((v, z), FSN_eq_with_stim)

# intersection with the dvf=0 fold curves: 
for v_value, z_value in zip(fold_lines[0], fold_lines[1]):
    if math.isclose(reduced_nullcline(v_value, z_value), 0, abs_tol = 0.0000001):
        FSN_values.append([v_value, z_value])

if FSN_values:
    FSN_z_value = [FSN_values[0][1]]
    FSN_v_value = [FSN_values[0][0]]

# reduced system: we have dvdtau = dyf*dydt + dzfdzdt, dzdtau = -dvfdzdt
# FSN_eq = FSN_eq.subs(y, phi)
# this setup is specifically used in the integrators.
FSN_eq_func_v = FSN_eq.subs(sympy_params_no_I)
FSN_eq_func_v = sp.lambdify((v, z, I), FSN_eq_func_v)
FSN_eq_z = -d_v_f*(dzdt)
FSN_eq_func_z = FSN_eq_z.subs(y, phi)
FSN_eq_func_z = FSN_eq_func_z.subs(sympy_params_no_I)
FSN_eq_func_z = sp.lambdify((v, z, I), FSN_eq_func_z)

# ODE problem setup, to be used in the integrator.
f_func = sp.lambdify((v, y, z, I), dvdt.subs(sympy_params_no_I))
h_func = sp.lambdify((v, y, I), dydt.subs(sympy_params_no_I))
g_func = sp.lambdify((v, z, I), dzdt.subs(sympy_params_no_I))

# setup for the original (non-dimensionalised) system
def nondim_ivp(t, X, I):
    v, y, z = X
    return[f_func(v, y, z, I), h_func(v, y, I), g_func(v, z, I)]

# setup for the layer problem
def nondim_layer_ivp(t, v, y, z, I):
    return[f_func(v, y, z, I)]

# setup for the reduced system
def reduced_system(X):
    v, z = X
    return[FSN_eq_func_v(v, z, I_stim), FSN_eq_func_z(v, z, I_stim)]

# setup for the reduced system integrator.
def reduced_system_ivp(t, X, stimulus):
    v, z = X
    dvdtau = FSN_eq_func_v(v, z, stimulus)
    dzdtau = FSN_eq_func_z(v, z, stimulus)
    return[dvdtau, dzdtau]

# To patch the GSPT flow together 
# TODO: work on this function tonight... it looks kinda yuck atm.
def patch_flow(eq, init_coords, plot_phi, projected_folds, given_folds, drop_up, drop_off, stimulus):
    t_span = [0, 10000]
    t_span_longer = [0, 500000]
    t_eval = np.linspace(t_span[0], t_span[1], points)
    t_eval_longer = np.linspace(t_span_longer[0], t_span_longer[1], points*1000)

    solution_v = []
    solution_z = []
    solution_y = []
    # see where the kick 'ends up' under the layer problem flow. 
    kick_up = solve_ivp(nondim_layer_ivp, t_span, [init_coords[0]], method=int_method, args=(init_coords[1], init_coords[2], stimulus), t_eval=t_eval)

    last_coords = [kick_up.y[0][-1], init_coords[1], init_coords[2]]
    solution_v.extend(kick_up.y[0])
    solution_z.extend([init_coords[2] for i in range(points)])
    solution_y.extend([init_coords[1] for i in range(points)])

    # Now we check, which 'part' of the critical manifold did the kicked point land on?
    # if it is less than the v value of the lower fold, then it must have landed on the 
    # lower portion of the manifold, and we only need to see what it does (likely follows to the equilibrium)
    if (kick_up.y[0][-1] <= given_folds[1]):
        # we can just compute how it follows the flow on the lower branch of the manifold.
        init_lower = [kick_up.y[0][-1], init_coords[2]]
        reduced_lower = solve_ivp(reduced_system_ivp, t_span, init_lower, args=(stimulus,), t_eval=t_eval, method=int_method)
        reduced_lower.y = reduced_lower.y[:, reduced_lower.y[0] <= given_folds[1]]

        solution_v.extend(reduced_lower.y[0])
        solution_z.extend(reduced_lower.y[1])
        solution_y.extend(plot_phi(reduced_lower.y[0], reduced_lower.y[1]))

        return [np.array(solution_v), np.array(solution_y), np.array(solution_z)]
                
    check = (not math.isclose(last_coords[0], eq)) or (not math.isclose(last_coords[1], plot_yinf(eq))) or (not math.isclose(last_coords[2], plot_zinf(eq)))
    i=0
    while check and i < 10:
        # if this is not the first-spike: we need to compute where the flow lands on the upper branch. 
        if not math.isclose(last_coords[0], kick_up.y[0][-1]):
            solution_v.extend(drop_up)
            solution_z.extend([reduced_lower.y[1][-1] for i in range(points)])
            solution_y.extend([plot_phi(reduced_lower.y[0][-1], reduced_lower.y[1][-1]) for i in range(points)])
            last_coords = [drop_up[-1], solution_y[-1], solution_z[-1]]

        # compute the reduced flow on the upper branch of the critical manifold
        reduced_upper = solve_ivp(reduced_system_ivp, t_span_longer, [last_coords[0], last_coords[2]], args=(stimulus,), t_eval=t_eval_longer, method=int_method)
        reduced_upper.y = reduced_upper.y[:, reduced_upper.y[0] >= given_folds[0]]

        solution_v.extend(reduced_upper.y[0])
        solution_z.extend(reduced_upper.y[1])
        solution_y.extend(plot_phi(reduced_upper.y[0], reduced_upper.y[1]))

        # Now the flow drops off the regular jump point.
        solution_v.extend(drop_off)
        solution_z.extend([reduced_upper.y[1][-1] for i in range(points)])
        solution_y.extend([plot_phi(reduced_upper.y[0][-1], reduced_upper.y[1][-1]) for i in range(points)])

        # from the drop off, compute the 
        # on the lower fold, compute the reduced flow again
        init_lower = [projected_folds[0], reduced_upper.y[1][-1]]
        reduced_lower = solve_ivp(reduced_system_ivp, t_span_longer, init_lower, args=(stimulus,) ,t_eval=t_eval_longer, method=int_method)
        reduced_lower.y = reduced_lower.y[:, reduced_lower.y[0] <= given_folds[1]]

        solution_v.extend(reduced_lower.y[0])
        solution_z.extend(reduced_lower.y[1])
        solution_y.extend(plot_phi(reduced_lower.y[0], reduced_lower.y[1]))
        
        last_coords = [solution_v[-1], solution_y[-1], solution_z[-1]]
        i=i+1
        check = (not math.isclose(last_coords[0], eq, abs_tol= 0.0001)) or (not math.isclose(last_coords[1], plot_yinf(eq), abs_tol= 0.0001)) or (not math.isclose(last_coords[2], plot_zinf(eq), abs_tol= 0.0001))

    return [np.array(solution_v), np.array(solution_y), np.array(solution_z)]

# TODO: fix this for the 'new' way I comput the folds
# the fold lines are computed with stimulus
def compute_projected_folds(fold_lines, stimulus, phi, plot_phi):
    projected_folds = []
    v_projected_fold = []
    z_projected_folds = fold_lines[1]
    phi_vals = plot_phi(fold_lines[0], fold_lines[1])
    # now we have z values, and y values we want the v values (I think)
    for z_val in z_projected_folds: 
        for val in phi_vals:
            phi_subs = phi.subs(z, z_val)
            # we want all v such that z and phi 'work' hmmmmm
            # i need to have a think through how i do this
    return projected_folds

def compute_invariant_manifolds(FSN_v_value, FSN_z_value):
    init = [float(FSN_v_value), float(FSN_z_value)]
    epsilon = 0.01
    t_inf = 0
    stable_inv_manifold = []
    unstable_inv_manifold = []

    Jac = jacobian(reduced_system, init)
    Jac = np.select([Jac.df > 1e-14, Jac.df < -1e-14], [Jac.df, Jac.df], 0)
    eigs, eigvecs = np.linalg.eig(Jac) 

    for i in range(2):
        if eigs[i] > 0:
            stable_eigvec = eigvecs[i]
        elif eigs[i] < 0:
            unstable_eigvec = eigvecs[i]

    t_backwards = [1800000, t_inf]
    t_eval_backwards = np.linspace(t_backwards[0], t_backwards[1], points*100)

    t_forwards = [t_inf, 1800000]
    t_eval_forwards = np.linspace(t_forwards[0], t_forwards[1], points*100)

    # Modify initial conditions to shoot in stable direction
    try:
        perturbed_stable = [[init[i] + 0.001 * stable_eigvec[i] for i in range(2)], [init[i] - 0.005 * stable_eigvec[i] for i in range(2)]]
        perturbed_unstable = [[init[i] + 0.05 * unstable_eigvec[i] for i in range(2)], [init[i] - 0.001 * unstable_eigvec[i] for i in range(2)]]
        stable_inv_manifold.append(solve_ivp(reduced_system_ivp, t_backwards, perturbed_stable[0], t_eval = t_eval_backwards, method=int_method, args=(I_stim,)))
        stable_inv_manifold.append(solve_ivp(reduced_system_ivp, t_backwards, perturbed_stable[1], t_eval = t_eval_backwards, method=int_method, args=(I_stim,)))
        unstable_inv_manifold.append(solve_ivp(reduced_system_ivp, t_forwards, perturbed_unstable[0], t_eval = t_eval_forwards, method=int_method, args=(I_stim,)))
        unstable_inv_manifold.append(solve_ivp(reduced_system_ivp, t_forwards, perturbed_unstable[1], t_eval = t_eval_forwards, method=int_method, args=(I_stim,)))
    except:
        print("something went wrong with the solver")
        stable_eig = 0
    
    return [stable_inv_manifold, unstable_inv_manifold]
    

class NonDimensionalThreeDim(System):
    def run():
        # attempting to make a dent at shooting in the direction of the stable invariant manifold - to get the 
        # (conjectured) separatrix that we will need to cross.
        if FSN_values:
            [stable_inv_manifold, unstable_inv_manifold] = compute_invariant_manifolds(FSN_v_value[0], FSN_z_value[0])

        z_vals = np.linspace(-0.3, 0.3, points) 
        v_vals = np.linspace(-0.9, 0.5, points)
        V, Z = np.meshgrid(v_vals, z_vals)

        # first: we plot the 'nullclines' with associated equilbria and folded saddle points
        plt.figure(figsize=(10, 6))
        for eq in equilibria:
            plt.plot(plot_zinf(eq), eq, marker='o', color='red', label="Equilibrium", zorder=10)
        plt.contour(Z, V, FSN_eq_func_v(V, Z, I_stim), levels=[0], linestyles = '--',colors='black', label = "v nullclines")
        plt.contour(Z, V, FSN_eq_func_z(V, Z, I_stim), levels=[0], color='black', label = "z nullclines")     
        if FSN_values:
            plt.plot(FSN_z_value[0], FSN_v_value[0], color='green', marker='o', label = "Folded Singularity")  
        plt.ylabel('v')
        plt.xlabel('z')
        plt.ylim(-0.9, 0.6)
        plt.xlim(-0.3,0.3)
        plt.title("Nullcline Plot with Associated Equilibria for De-singularised Reduced Problem")
        plt.grid()
        plt.legend()
        plt.show()

        # next we plot all the invariant manifolds for the folded singularity, with the folds and Folded Saddle included
        # in this case, I colour the alternate (in-between the folds) the stability of the manifold. 
        if stable_inv_manifold:
            plt.plot((stable_inv_manifold[0]).y[1], (stable_inv_manifold[0]).y[0], color='blue', label= "Stable invariant curve")
            plt.plot((stable_inv_manifold[1]).y[1], (stable_inv_manifold[1]).y[0], '--', color='red', label= "Flipped stable invariant curve") 
        if unstable_inv_manifold:
            plt.plot((unstable_inv_manifold[0]).y[1], (unstable_inv_manifold[0]).y[0], color='red', label="Unstable invariant curve")
            plt.plot((unstable_inv_manifold[1]).y[1], (unstable_inv_manifold[1]).y[0], '--', color='blue', label="Flipper unstable invariant curve")
        if FSN_values:
            plt.plot(FSN_z_value[0], FSN_v_value[0], color='green', marker='o', label = "Folded Singularity", zorder=10)
        plt.plot(fold_lines[1], fold_lines[0], color='black', label="Lower Fold line")
        plt.ylabel('v')
        plt.xlabel('z')
        plt.ylim(-1, 0.05)
        plt.xlim(-0.3,0.3)
        plt.title("Invariant Manifolds of Folded Singularity")
        plt.grid()
        if FSN_z_value and FSN_z_value[0] > 0:
            plt.legend(loc='lower left')
        else:
            plt.legend(loc='lower right')
        plt.show()

        # obtain the kicked equilibria starting when the stimulus is 0.
        initial_conditions = [og_equilibria[0], float(plot_yinf_og(og_equilibria[0])), float(plot_zinf_og(og_equilibria[0]))]

        # projected fold line to the opposite branch
        # projected_folds = compute_projected_folds(folds, I_stim)

        # get all my plot features together, making sure to only plot when z, y are between 0 and 1.
        v_vals, z_vals = np.meshgrid(v_vals, z_vals)
        y_vals = plot_phi(v_vals, z_vals)
        y_vals = np.where((y_vals >= 0) & (y_vals <= 0.3), y_vals, np.nan)
        fold_y_vals = plot_phi(np.asarray(fold_lines[0]), np.asarray(fold_lines[1]))
        fold_y_vals = np.where((fold_y_vals >= 0) & (fold_y_vals <= 0.3), fold_y_vals, np.nan)

        # folds along the manifold
        # y_upper_fold = plot_phi(folds[0], z_vals)
        # y_upper_fold = np.where((y_upper_fold >= 0) & (y_upper_fold <= 0.3), y_upper_fold, np.nan)
        # y_lower_fold = plot_phi(folds[1], z_vals)
        # y_lower_fold = np.where((y_lower_fold >= 0) & (y_lower_fold <= 0.3), y_lower_fold, np.nan)

        # drop_off = np.linspace(folds[0], projected_folds[0], points)
        # drop_up = np.linspace(folds[1], projected_folds[1], points)
        
        # stable_inv_manifold.y = stable_inv_manifold.y[:, stable_inv_manifold.y[0] > projected_folds[0]]
        
        # transient_sols = patch_flow(equilibria[0], initial_conditions, plot_phi, projected_folds, folds, drop_up, drop_off, I_stim)
        t_span = [0, 20000]
        t_eval = np.linspace(t_span[0], t_span[1], points*1000)
        actual_sol = solve_ivp(nondim_ivp, t_span, initial_conditions, method=int_method, t_eval=t_eval, args=(I_stim,))

        plt.plot(actual_sol.t, actual_sol.y[0], color='black', linewidth=1.5, label="transient spike time trace")     
        plt.xlabel('t (non-dimensional time)')
        plt.ylabel('v (non-dimensional voltage term)')
        plt.title("Transient Spiking Behaviour in Type 3 Neuron")
        plt.grid()
        plt.legend()
        plt.show()

        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')

        # Now we can plot the slow manifold (y = phi(v, z)), with equilibria and fold lines included.
        ax.plot_surface(y_vals, z_vals, v_vals, label='Critical manifold', color='sandybrown', alpha=0.9)
        for eq in equilibria:
            ax.scatter(plot_yinf(eq), plot_zinf(eq), eq, marker='o', color='red', label='equilibrium')  
        if FSN_values and FSN_z_value[0] >= 0:
            ax.scatter(float(plot_phi(FSN_v_value[0], FSN_z_value[0])), float(FSN_z_value[0]), float(FSN_v_value[0]), marker='o', color='green', label='FSN')
            # ax.plot(plot_phi(stable_inv_manifold.y[0], stable_inv_manifold.y[1]), stable_inv_manifold.y[1], stable_inv_manifold.y[0], color='blue', label= "Stable invariant curve") 
        # plot the fold lines, and the projection of the fold onto the opposite 'branch'.
        plt.plot(fold_y_vals, np.asarray(fold_lines[1]), np.asarray(fold_lines[0]), color='black', ls='--', linewidth= 2, label='folds?')
        # ax.plot(y_upper_fold, z_vals, folds[0], color='black', ls = '--', linewidth=2, label="Upper Fold")
        # ax.plot(y_lower_fold, z_vals, folds[1], color='black', ls = '--', linewidth =2, label= "lower Fold") 
        # ax.plot(y_upper_fold, z_vals, projected_folds[0], color='darkblue', ls = '--', linewidth =2, label= "Projected Upper fold")
        # ax.plot(y_lower_fold, z_vals, projected_folds[1], color='darkblue', ls = '--', linewidth =2, label= "Projected Lower fold")                        
        # ax.plot(transient_sols[1], transient_sols[2], transient_sols[0], color='darkmagenta', linewidth=1.5, alpha=1, label="Solution")  
        # ax.scatter(initial_conditions[1], initial_conditions[2], initial_conditions[0], marker='o', color='pink', label='kicked point') 
        # ax.plot(actual_sol.y[1], actual_sol.y[2], actual_sol.y[0], color='black', lw=2)     
        ax.set_xlabel('y')
        ax.set_ylabel('z')
        ax.set_zlabel('v')
        ax.set_xlim(0,0.3)
        ax.set_ylim(0,0.3)
        ax.set_zlim(-1, 0.5)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True)
        mplstyle.use('fast')
        #plt.savefig("critical manifold full plot")
        plt.show()

    def move():
        Z_vals = np.linspace(0, 0.6, points) 
        V_vals = np.linspace(-0.9, 0.5, points)
        v_vals, z_vals = np.meshgrid(V_vals, Z_vals)

        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')

        plt.subplots_adjust(bottom=0.25) 

        ax.set_xlim(0, 0.4)
        ax.set_ylim(0, 0.4)
        ax.set_zlim(-1, 0.5)        
        ax.set_xlabel('y')
        ax.set_ylabel('z')
        ax.set_label('v')
        ax.set_title('Changes to spiking behaviour as stimulus is varied')

        initial_conditions = [og_equilibria[0], plot_yinf(og_equilibria[0]), plot_zinf(og_equilibria[0])]

        y_vals = plot_phi(v_vals, z_vals)
        y_vals = np.where((y_vals >= 0) & (y_vals <= 0.3), y_vals, np.nan)

        y_upper_fold = plot_phi(folds[0], Z_vals)
        y_upper_fold = np.where((y_upper_fold >= 0) & (y_upper_fold <= 0.3), y_upper_fold, np.nan)
        y_lower_fold = plot_phi(folds[1], Z_vals)
        y_lower_fold = np.where((y_lower_fold >= 0) & (y_lower_fold <= 0.3), y_lower_fold, np.nan)

        projected_folds = compute_projected_folds(folds, I_stim)
        drop_off = np.linspace(folds[0], projected_folds[0], points)
        drop_up = np.linspace(folds[1], projected_folds[1], points)

        transient_sols = patch_flow(equilibria[0], initial_conditions, plot_phi, projected_folds, folds, drop_up, drop_off, I_stim)

        ax.plot_surface(y_vals, z_vals, v_vals, label='Critical manifold', color='sandybrown', alpha=0.8)
        ax.scatter(plot_yinf(equilibria[0]), plot_zinf(equilibria[0]), equilibria[0], color='red', label='Equilibria')
        ax.scatter(plot_phi(folds[1], FSN_z_value[0]), FSN_z_value[0], folds[1], color='green', label='Folded Sigularity')
        ax.scatter(initial_conditions[1], initial_conditions[2], initial_conditions[0], color='pink', label="Initial Condition")
        ax.plot(y_lower_fold, Z_vals, folds[1], color='black', label='Lower Fold', lw=1, ls='--')
        ax.plot(y_upper_fold, Z_vals, folds[0], color='black', label='Upper Fold', lw=1, ls='--')
        ax.plot(y_lower_fold, Z_vals, projected_folds[1], color='darkblue', ls = '--', linewidth =2, label= "Projected Lower fold")   
        ax.plot(y_upper_fold, Z_vals, projected_folds[0], color='darkblue', ls = '--', linewidth =2, label= "Projected Lower fold")                                             
        ax.plot(transient_sols[1], transient_sols[2], transient_sols[0], color='darkmagenta', linewidth=1.5, alpha=1, label="Solution")  

        I_text = ax.text(0, 0, 0, '', transform=ax.transAxes,
                        verticalalignment='top', fontsize=12,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # slider:
        ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03])
        slider = Slider(ax_slider, 'I', valmin=0.0, valmax=100/(kv*20), valinit=I_stim, valstep=0.0001, facecolor="black", edgecolor="black")  

        slider.label.set_fontsize(12)
        slider.valtext.set_fontweight("bold")

        def update(val):
            I_current = slider.val

            ax.clear()
            ax.scatter(plot_yinf(og_equilibria[0]), plot_zinf(og_equilibria[0]), og_equilibria[0], color='pink', label="Initial Condition")

            current_phi = phi_subs_no_I.subs(I, I_current)
            plot_current_phi = sp.lambdify((v,z), current_phi)
            current_equilibria = util.find_equilibria(eq_curve_subs, v, eq_guesses, I_current)
            current_folds = util.find_equilibria(eq_fold_subs, v, fold_guesses, I_current)
            current_fold_proj = compute_projected_folds(current_folds, I_current)
            current_fall = np.linspace(current_folds[0], current_fold_proj[0], points)
            current_jump = np.linspace(current_folds[1], current_fold_proj[1], points)
            transient_sols = patch_flow(current_equilibria[0], initial_conditions, plot_current_phi, current_fold_proj, current_folds, current_jump, current_fall, I_current)


            FSN_val_changing = FSN_eq_subs.subs(v, current_folds[1]) if current_folds else None
            if FSN_val_changing is not None:
                FSN_val_changing = FSN_val_changing.subs(I, I_current)
                current_FSN = sp.solve(FSN_val_changing, z)
            else:
                current_FSN = None
            
            y_new = plot_current_phi(v_vals, z_vals)
            y_new = np.where((y_new >= 0) & (y_new <= 0.3), y_new, np.nan)
            
            ax.plot_surface(y_new, z_vals, v_vals, label='Critical manifold', color='sandybrown', alpha=0.8)

            # Update Equilibria
            if current_equilibria:
                y_new = np.atleast_1d(plot_yinf(current_equilibria[0]))
                z_new = np.atleast_1d(plot_zinf(current_equilibria[0]))
                v_new = np.atleast_1d(current_equilibria[0])
                ax.scatter(y_new, z_new, v_new, color='red', label='Equilibria')
            else:
                pass

            # Update the fold lines, and projected fold lines
            if current_folds:                
                y_upper_fold = plot_current_phi(current_folds[0], Z_vals)
                y_upper_fold = np.where((y_upper_fold >= 0) & (y_upper_fold <= 0.3), y_upper_fold, np.nan)
                y_lower_fold = plot_current_phi(current_folds[1], Z_vals)
                y_lower_fold = np.where((y_lower_fold >= 0) & (y_lower_fold <= 0.3), y_lower_fold, np.nan)

                ax.plot(y_lower_fold, Z_vals, current_folds[1], color='black', label='Lower Fold', lw=2, ls='--')
                ax.plot(y_upper_fold, Z_vals, current_folds[0], color='black', label='Upper Fold', lw=2, ls='--')
                ax.plot(y_lower_fold, Z_vals, current_fold_proj[1], color='darkblue', ls = '--', linewidth =2, label= "Projected Lower fold")   
                ax.plot(y_upper_fold, Z_vals, current_fold_proj[0], color='darkblue', ls = '--', linewidth =2, label= "Projected Upper fold")
                ax.plot(transient_sols[1], transient_sols[2], transient_sols[0], color='darkmagenta', linewidth=1.5, alpha=1, label="Solution")
            else:
                pass

            # update the FSN
            if current_FSN:
                y_val = plot_current_phi(current_folds[1], current_FSN[0])
                y_new = np.atleast_1d(np.float64(y_val))
                z_new = np.atleast_1d(np.float64(current_FSN[0]))
                v_new = np.atleast_1d(current_folds[1])
                ax.scatter(y_new, z_new, v_new, color='green', label='Folded Sigularity')
            else:
                pass

            # --- Text update ---
            I_text.set_text(f'Original system I = {I_current*kv*20:.1f}')
            mplstyle.use('fast')

            fig.canvas.draw_idle()

        slider.on_changed(update)
        ax.legend()
        plt.show()














    def move_old():
        # Now I want to animate how the equilibria move.
        eq_curve_no_I = eq_curve.subs(sympy_params_no_I)
        fold_curve_no_I = eq_fold.subs(sympy_params_no_I)
        FSN_eq_no_I = FSN_eq.subs(sympy_params_no_I)

        fig, ax = plt.subplots()

        ax.set_xlim(-0.5, 0.5)     
        ax.set_ylim(-1, 0.5)  
        ax.set_xlabel('z')
        ax.set_ylabel('v')
        ax.set_title('Movement of Folds & Equilibria as Stimulus is Varied')
        equilibria_scatter = ax.scatter([], [], c='red', s=15, zorder=5, label='Equilibria')
        FSN_scatter = ax.scatter([], [], c='green', s=15, zorder=5, label='FSN')
        plot_lower_fold, = ax.plot([],[], color='black', label='fold', lw=0.5)
        plot_upper_fold, = ax.plot([],[], color='black', label='fold', lw=0.5)
        I_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, 
                verticalalignment='top', fontsize=12,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        def update_data(frame):
            I_current = I_stim + (frame*0.1)/(kv*20)
            current_equilibria = util.find_equilibria(eq_curve_no_I, v, eq_guesses, I_current)
            current_folds = util.find_equilibria(fold_curve_no_I, v, fold_guesses, I_current)
            FSN_val_changing = FSN_eq_no_I.subs(v, current_folds[1])
            FSN_val_changing = FSN_val_changing.subs(I, I_current)
            current_FSN = sp.solve(FSN_val_changing, z)

            #if we find equilibria, update the graph.
            if current_equilibria:
                # Calculate z values for each equilibrium
                current_z_vals = []
                for eq in current_equilibria:
                    current_z_vals.append(float(plot_zinf(eq)))
                    # Update scatter plot with new positions
                plot_data = list(zip(current_z_vals, current_equilibria))
                equilibria_scatter.set_offsets(plot_data)
            else:
                # No equilibria found
                equilibria_scatter.set_offsets([])
            
            if current_folds:
                current_upper_fold = np.full(points, current_folds[0])
                current_lower_fold = np.full(points, current_folds[1])
                plot_upper_fold.set_data(z_vals, current_upper_fold)
                plot_lower_fold.set_data(z_vals, current_lower_fold)

            else:
                plot_upper_fold.set_data([],[])
                plot_lower_fold.set_data([],[])

            #if we find a FSN, update the graph.
            if current_FSN:
                # the z value is the given value when solving, ans v value is the fold value. 
                plot_data_FSN = [(current_FSN[0], current_folds[1])]
                FSN_scatter.set_offsets(plot_data_FSN)
            else:
                # No equilibria found
                FSN_scatter.set_offsets([])
            
            I_text.set_text(f'I = {I_current:.2f}')

            return equilibria_scatter, plot_lower_fold, plot_upper_fold, I_text


        anim = animation.FuncAnimation(fig = fig,
                                       func = update_data,
                                       frames=1000,
                                       interval = 30)

        plt.show()

        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)  # make space for slider

        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-1, 0.5)
        ax.set_xlabel('z')
        ax.set_ylabel('v')
        ax.set_title('Movement of Folds & Equilibria as Stimulus is Varied')
        ax.plot(np.zeros(points), np.linspace(-1, 0.5, points), color ='blue', lw=0.5)

        equilibria_scatter = ax.scatter([], [], c='red', s=15, zorder=5, label='Equilibria')
        FSN_scatter = ax.scatter([], [], c='green', s=15, zorder=5, label='Folded Sigularity')
        plot_lower_fold, = ax.plot([],[], color='black', label='Lower Fold', lw=1)
        plot_upper_fold, = ax.plot([],[], color='black', label='Upper Fold', lw=1)
        I_text = ax.text(0.02, 0.98, '', transform=ax.transAxes,
                        verticalalignment='top', fontsize=12,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # alternate option for a slider:
        ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03])
        slider = Slider(ax_slider, 'I', valmin=0.0, valmax=100/(kv*20), valinit=0.5, valstep=0.0001, facecolor="black", edgecolor="black")  

        slider.label.set_fontsize(12)
        slider.valtext.set_fontweight("bold")

        # --- Update function ---
        def update(val):
            I_current = slider.val

            current_equilibria = util.find_equilibria(eq_curve_no_I, v, eq_guesses, I_current)
            current_folds = util.find_equilibria(fold_curve_no_I, v, fold_guesses, I_current)

            FSN_val_changing = FSN_eq_no_I.subs(v, current_folds[1]) if current_folds else None
            if FSN_val_changing is not None:
                FSN_val_changing = FSN_val_changing.subs(I, I_current)
                current_FSN = sp.solve(FSN_val_changing, z)
            else:
                current_FSN = None

            # --- Equilibria update ---
            if current_equilibria:
                current_z_vals = [float(plot_zinf(eq)) for eq in current_equilibria]
                plot_data = list(zip(current_z_vals, current_equilibria))
                equilibria_scatter.set_offsets(plot_data)
            else:
                equilibria_scatter.set_offsets([])

            # --- Folds update ---
            if current_folds:
                current_upper_fold = np.full(points, current_folds[0])
                current_lower_fold = np.full(points, current_folds[1])
                plot_upper_fold.set_data(z_vals, current_upper_fold)
                plot_lower_fold.set_data(z_vals, current_lower_fold)
            else:
                plot_upper_fold.set_data([],[])
                plot_lower_fold.set_data([],[])

            # --- FSN update ---
            if current_FSN:
                plot_data_FSN = [(current_FSN[0], current_folds[1])]
                FSN_scatter.set_offsets(plot_data_FSN)
            else:
                FSN_scatter.set_offsets([])

            # --- Text update ---
            I_text.set_text(f'Original system I = {I_current*kv*20:.1f}')

            fig.canvas.draw_idle()  # redraw

            slider.on_changed(update)
            ax.legend()
            plt.show()
