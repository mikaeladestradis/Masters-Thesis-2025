import sympy as sp
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d.axes3d import Axes3D
import itertools
import util
from scipy.integrate import solve_ivp
from scipy.differentiate import jacobian
from systems.system import System
import matplotlib.style as mplstyle
import three_dimensional_system_pars
from matplotlib import animation
import warnings

warnings.filterwarnings("ignore")
points = 1000
kv = 100
int_method = "LSODA"

#for presentation - no transient (~20), one transient (~50), two transient (~80) three? #82.35 periodic 83
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

# parameter definition: a full list of parameters is given in the three_dimensional_system_pars file.
params = three_dimensional_system_pars.way_more_than_three
I_stim = params.get("I")

# different parameter setups. 
sympy_params = {sp.symbols(k): v for k, v in params.items()}
sympy_params_no_I = {sp.symbols(k): v for k, v in params.items() if k != 'I'}

# this is to compute and plot the function y = phi
phi_subs = phi.subs(sympy_params)
phi_subs_no_I = phi.subs(sympy_params_no_I)
plot_phi = sp.lambdify((v, z), phi_subs)
phi_subs_zero_I = phi_subs_no_I.subs("I", 0)
plot_phi_I_zero = sp.lambdify((v, z), phi_subs_zero_I)

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
eq_fold_no_z = eq_fold.subs(z, zinf)
eq_fold_no_z = eq_fold_no_z.subs(sympy_params_no_I)
fold_lines = util.find_fold_lines(eq_fold, v_vals, z_vals, [v, z], I_stim, sympy_params_no_I)
fold_v_vals = util.find_equilibria(eq_fold_no_z, v, fold_guesses, I_stim)

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

ramp = np.linspace(0, (100/(20*100)) , 29001)

def step_protocol(time, injected_current):
    if time <= 1000:
        return 0
    else:
        return injected_current

# setup for the original (non-dimensionalised) system
def nondim_ivp(t, X, injected_current):
    I = step_protocol(t, injected_current)

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

# TODO: fix this for the 'new' way I compute the folds
def compute_projected_folds(fold_lines, v_vals): # note for the time being this should work, but it takes a hella long time. 
    projected_folds_v = [] 
    phi_vals = [] 
    for i in range(len(fold_lines[0])): 
        phi_vals.append(plot_phi(fold_lines[0][i], fold_lines[1][i]))

    for (z_index, z_val), v_val in itertools.product(enumerate(fold_lines[1]), v_vals): 
        for phi_val in phi_vals: 
            print(phi_val)
            if math.isclose(plot_phi(v_val, z_val), phi_val, abs_tol=0.0001): 
                if (fold_lines[0][z_index] != v_val): 
                    projected_folds_v.append(v_val) 
    
    return projected_folds_v


def compute_invariant_manifolds(FSN_v_value, FSN_z_value):
    init = [float(FSN_v_value), float(FSN_z_value)]
    epsilon = 0.001
    t_inf = 0
    stable_inv_manifold = []
    unstable_inv_manifold = []

    Jac = jacobian(reduced_system, init)
    Jac = np.select([Jac.df > 1e-14, Jac.df < -1e-14], [Jac.df, Jac.df], 0)
    eigs, eigvecs = np.linalg.eig(Jac) 
    print(eigs)

    for i in range(2):
        if eigs[i] > 0:
            stable_eigvec = eigvecs[i]
        elif eigs[i] < 0:
            unstable_eigvec = eigvecs[i]

    t_backwards = [1800000, t_inf]
    t_eval_backwards = np.linspace(t_backwards[0], t_backwards[1], points*10000)

    t_forwards = [t_inf, 1800000]
    t_eval_forwards = np.linspace(t_forwards[0], t_forwards[1], points*1000)

    # Modify initial conditions to shoot in stable direction
    try:
        perturbed_stable = [[init[i] + 0.00001 * stable_eigvec[i] for i in range(2)], [init[i] - 0.0015 * stable_eigvec[i] for i in range(2)]]
        perturbed_unstable = [[init[i] + 0.001 * unstable_eigvec[i] for i in range(2)], [init[i] - 0.001 * unstable_eigvec[i] for i in range(2)]]
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
        v_vals = np.linspace(-0.9, 0.6, points)
        V, Z = np.meshgrid(v_vals, z_vals)

        # first: we plot the 'nullclines' with associated equilbria and folded saddle points
        plt.figure(figsize=(10, 6))
        for eq in equilibria:
            plt.plot(plot_zinf(eq), eq, marker='o', color='red', label="True Equilibrium", zorder=10)
        plt.contour(Z, V, FSN_eq_func_v(V, Z, I_stim), levels=[0], linestyles = '--',colors='black', label = "v nullclines")
        plt.plot([], [], 'k--', label='v nullcline')
        plt.plot([], [], 'k-', label='w nullcline')
        plt.contour(Z, V, FSN_eq_func_z(V, Z, I_stim), levels=[0], color='black', label = "z nullclines")     
        if FSN_values:
            plt.plot(FSN_z_value[0], FSN_v_value[0], color='green', marker='o', label = "Folded Singularity", zorder=10)  
        plt.ylabel('v')
        plt.xlabel('z')
        plt.ylim(-0.9, 0.6)
        plt.xlim(-0.3,0.3)
        plt.title("Nullcline Plot for De-singularised Reduced Problem")
        plt.grid()
        plt.legend()
        plt.show()

        # next we plot all the invariant manifolds for the folded singularity, with the folds and Folded Saddle included
        # in this case, I colour the alternate (in-between the folds) the stability of the manifold. 
        if FSN_values:
            if stable_inv_manifold:
                plt.plot((stable_inv_manifold[0]).y[1], (stable_inv_manifold[0]).y[0], color='purple', label= "canard solution")
                plt.plot((stable_inv_manifold[1]).y[1], (stable_inv_manifold[1]).y[0], color='purple') 
            if unstable_inv_manifold:
                plt.plot((unstable_inv_manifold[0]).y[1], (unstable_inv_manifold[0]).y[0], color='green', label="faux canard solution")
                plt.plot((unstable_inv_manifold[1]).y[1], (unstable_inv_manifold[1]).y[0], color='green')
            # if FSN_values:
            #     plt.plot(FSN_z_value[0], FSN_v_value[0], color='green', marker='o', label = "Folded Singularity", zorder=10)
            for eq in equilibria:
                plt.plot(plot_zinf(eq), eq, marker='o', color='red', label="True Equilibrium", zorder=10)
            plt.plot(fold_lines[1], fold_lines[0], color='black', label=r"$F^-$")
            plt.ylabel('v')
            plt.xlabel('z')
            plt.ylim(-1, 0.05)
            plt.xlim(-0.3,0.3)
            plt.title("Invariant Manifolds of Folded Singularity")
            plt.grid()
            if FSN_z_value and FSN_z_value[0] > 0:
                plt.legend(loc='lower left', bbox_to_anchor=(0, 0))
            else:
                plt.legend(loc='lower right', bbox_to_anchor=(1, 0))
            plt.show()

        # obtain the kicked equilibria starting when the stimulus is 0.
        initial_conditions = [og_equilibria[0], float(plot_yinf_og(og_equilibria[0])), float(plot_zinf_og(og_equilibria[0]))]

        # get all my plot features together, making sure to only plot when z, y are between 0 and 1.
        v_vals = np.linspace(-0.9, 0.6, points)
        z_vals = np.linspace(0, 0.3, points)
        v_vals, z_vals = np.meshgrid(v_vals, z_vals)

        # next we want to restrict to relevant domain (ie) y and z between 0 and 1)
        y_vals = plot_phi(v_vals, z_vals)
        y_vals = np.where((y_vals >= 0) & (y_vals <= 0.3), y_vals, np.nan)
        fold_y_vals = plot_phi(np.asarray(fold_lines[0]), np.asarray(fold_lines[1]))
        fold_y_vals = np.where((fold_y_vals >= 0) & (fold_y_vals <= 0.3), fold_y_vals, np.nan)
        fold_lines[1] = [val if ((0 <= val) & (val <= 0.3)) else np.nan for val in fold_lines[1]]
        if FSN_values and stable_inv_manifold:
            y_separatrix_values = plot_phi((stable_inv_manifold[0]).y[0], (stable_inv_manifold[0]).y[1])
            y_separatrix_values = np.where((y_separatrix_values >= 0) & (y_separatrix_values <= 0.3), y_separatrix_values, np.nan)
        
        # transient_sols = patch_flow(equilibria[0], initial_conditions, plot_phi, projected_folds, folds, drop_up, drop_off, I_stim)
        t_span = [0, 30000]
        t_eval = np.linspace(t_span[0], t_span[1], points*100)

        actual_sol = solve_ivp(nondim_ivp, t_span, initial_conditions, method=int_method, t_eval=t_eval, args=(I_stim,))

        step_protocol = [I_stim if val >= 1000 else 0 for val in t_eval]

        actual_sol = solve_ivp(
            nondim_ivp,
            t_span,
            initial_conditions,
            method=int_method,
            t_eval=t_eval,
            args=(I_stim,)
        )

        t = actual_sol.t
        v = actual_sol.y[0]

        frame_indices = np.arange(0, len(t), 10000)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,6), sharex=True)

        ax1.set_xlim(t_span)
        ax1.set_ylim(-0.9, 0.5)
        ax1.set_ylabel("v(t)")
        ax1.set_title("Transient Behaviour")

        voltage_line, = ax1.plot([], [], color="black", linewidth=1.5)

        ax2.set_xlim(t_span)
        ax2.set_ylim(0, -0.08)
        ax2.set_xlabel("t")
        ax2.set_ylabel("Injected Current")

        stim_line, = ax2.plot([], [], color="black", linewidth=1.5)

        def init():
            voltage_line.set_data([], [])
            stim_line.set_data([], [])
            return voltage_line, stim_line

        def update(frame):
            idx = frame_indices[frame]
            voltage_line.set_data(t[:idx], v[:idx])
            stim_line.set_data(t_eval[:idx], step_protocol[:idx])
            return voltage_line, stim_line

        anim = animation.FuncAnimation(
            fig,
            update,
            frames=len(frame_indices),
            init_func=init,
            interval=10,
            blit=True
        )

        # writer = animation.PillowWriter(fps=50)
        # anim.save("one-spike.gif", writer=writer)
        plt.show()

        fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

        axs[0].plot(actual_sol.t, actual_sol.y[0], color='black')
        axs[0].set_ylabel("v")
        axs[0].set_ylim(-0.9, 0.5)
        axs[0].set_title("Time trace and injected current")

        axs[1].plot(actual_sol.t, step_protocol, color='black')
        axs[1].set_xlabel("T")
        axs[1].set_ylim(0, 0.05)
        axs[1].set_ylabel(r"$I$")

        plt.tight_layout()
        plt.savefig("time-trace-three-spikes")
        plt.close()

        # fig = plt.figure(figsize=(10, 8), constrained_layout=True)
        # ax = fig.add_subplot(111, projection='3d')

        # upper_fold_v_val = fold_v_vals[0]
        # lower_fold_v_val = fold_v_vals[1]

        # upper_fold_v_vals = np.full(points, upper_fold_v_val)
        # lower_fold_v_vals = np.full(points, lower_fold_v_val)

        # upper_fold_y_vals = plot_phi(upper_fold_v_val, z_vals)
        # lower_fold_y_vals = plot_phi(lower_fold_v_val, z_vals)
        # upper_fold_y_vals = np.where((upper_fold_y_vals >= 0) & (upper_fold_y_vals <= 0.3), upper_fold_y_vals, np.nan)
        # lower_fold_y_vals = np.where((lower_fold_y_vals >= 0) & (lower_fold_y_vals <= 0.3), lower_fold_y_vals, np.nan)

        # upper_stable_branch = []
        # lower_stable_branch = []
        # unstable_branch = []
        # v_values = np.linspace(-0.9, 0.4, points)

        # for value in v_values:
        #     if (value <= lower_fold_v_val):
        #         lower_stable_branch.append(value)
        #     elif (value >= upper_fold_v_val):
        #         upper_stable_branch.append(value)
        #     else:
        #         unstable_branch.append(value)
        
        # z_values_upper = np.linspace(0, 0.3, len(upper_stable_branch))
        # z_values_lower = np.linspace(0, 0.3, len(lower_stable_branch))
        # z_values_middle = np.linspace(0, 0.3, len(unstable_branch))

        # V_VALS_UPPER, Z_VALS_UPPER = np.meshgrid(np.array(upper_stable_branch), z_values_upper)
        # V_VALS_LOWER, Z_VALS_LOWER = np.meshgrid(np.array(lower_stable_branch), z_values_lower)
        # V_VALS_MIDDLE, Z_VALS_MIDDLE = np.meshgrid(np.array(unstable_branch), z_values_middle)

        # upper_stable_phi = plot_phi(V_VALS_UPPER, Z_VALS_UPPER)
        # lower_stable_phi = plot_phi(V_VALS_LOWER, Z_VALS_LOWER)
        # unstable_phi = plot_phi(V_VALS_MIDDLE, Z_VALS_MIDDLE)

        # upper_stable_phi = np.where((upper_stable_phi >= 0) & (upper_stable_phi <= 0.3), upper_stable_phi, np.nan)
        # lower_stable_phi = np.where((lower_stable_phi >= 0) & (lower_stable_phi <= 0.3), lower_stable_phi, np.nan)
        # unstable_phi = np.where((unstable_phi >= 0) & (unstable_phi <= 0.3), unstable_phi, np.nan)

        # # Now we can plot the slow manifold (y = phi(v, z)), with equilibria and fold lines included.
        # # ax.plot_surface(upper_stable_phi, Z_VALS_UPPER, V_VALS_UPPER, color='blue', alpha=0.4, label=r"$S^{\pm}$")
        # # ax.plot_surface(lower_stable_phi, Z_VALS_LOWER, V_VALS_LOWER, color='blue', alpha=0.4)        
        # # ax.plot_surface(unstable_phi, Z_VALS_MIDDLE, V_VALS_MIDDLE, color='red', label=r"$U$", alpha=0.4)
        # # ax.plot_surface(np.full(len(upper_stable_branch), -0.05), Z_VALS_UPPER, V_VALS_UPPER, color='blue', alpha=0.4)
        # # ax.plot_surface(np.full(len(lower_stable_branch), -0.05), Z_VALS_LOWER, V_VALS_LOWER, color='blue', alpha=0.4)        
        # # ax.plot_surface(np.full(len(unstable_branch), -0.05), Z_VALS_MIDDLE, V_VALS_MIDDLE, color='red', alpha=0.4)
        # ax.plot_surface(y_vals, z_vals, v_vals, label=r'$M$', color='grey', alpha=0.4)
        # for eq in equilibria:
        #     ax.scatter(plot_yinf(eq), plot_zinf(eq), eq, marker='o', color='red', label='Equilibrium')  
        # # if FSN_values and FSN_z_value[0] >= 0:
        # #     ax.scatter(float(plot_phi(FSN_v_value[0], FSN_z_value[0])), float(FSN_z_value[0]), float(FSN_v_value[0]), marker='o', color='green', label='FSN')
        # # plot the fold lines, and the projection of the fold onto the opposite 'branch'.
        # plt.plot(upper_fold_y_vals, z_vals, upper_fold_v_vals, color='purple', ls='--', linewidth= 3, label=r"$F^+$")
        # plt.plot(lower_fold_y_vals, z_vals, lower_fold_v_vals, color='green', ls='--', linewidth= 3, label=r"$F^-$")
        # # plt.plot(np.full(points, -0.05), z_vals, upper_fold_v_vals, color='purple', ls='--', linewidth= 3)
        # # plt.plot(np.full(points, -0.05), z_vals, lower_fold_v_vals, color='green', ls='--', linewidth= 3)
        # ax.scatter(initial_conditions[1], initial_conditions[2], initial_conditions[0], marker='o', color='pink', label='Initial Condition') 
        # ax.plot(actual_sol.y[1], actual_sol.y[2], actual_sol.y[0], color='black', lw=2, label="Solution Trajectory")     
        # if FSN_values and stable_inv_manifold:
        #     ax.plot(y_separatrix_values, (stable_inv_manifold[0]).y[1], (stable_inv_manifold[0]).y[0], color='purple', label= "Canard Solution")
        # ax.set_xlabel('y')
        # ax.set_ylabel('z')
        # ax.set_zlabel('v')
        # ax.set_xlim(-0.05,0.3)
        # ax.set_ylim(0,0.3)
        # ax.set_zlim(-1, 0.5)
        # # ax.legend(loc='upper right', fancybox=True)
        # ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=True)
        # mplstyle.use('fast')
        # #plt.savefig("critical manifold full plot")
        # plt.show()
        
        # hold_z_init = plot_zinf(initial_conditions[0])
        # v_vals = np.linspace(-0.9, 0.5, 1000)
        # og_y_vals = plot_phi_I_zero(v_vals, hold_z_init)

        # plt.plot(og_y_vals, v_vals, color='grey')
        # plt.scatter(plot_yinf(initial_conditions[0]), initial_conditions[0], color='red', zorder=10, s=20)
        # plt.xlim(-0.1, 0.3)
        # plt.ylim(-0.9, 0.4)
        # plt.show()

        # hold_z_eq = plot_zinf(equilibria[0])
        # v_vals = np.linspace(-0.9, 0.5, 1000)
        # y_vals_2d = plot_phi(v_vals, hold_z_eq)

        # plt.plot(y_vals_2d, v_vals, color='grey')
        # plt.scatter(initial_conditions[1], initial_conditions[0], color='pink', zorder=10, s=20)
        # plt.scatter(plot_yinf(equilibria[0]), equilibria[0], color='red', zorder=10, s=20)
        # plt.plot(actual_sol.y[1], actual_sol.y[0], color='black')
        # plt.xlim(-0.05, 0.25)
        # plt.ylim(-0.9, 0.4)
        # plt.xlabel("y")
        # plt.ylabel('v')
        # plt.show()

        # fig, ax1 = plt.subplots()
        # ax1.scatter(initial_conditions[1], initial_conditions[0], color='pink', zorder=10, s=20)
        # ax1.set_xlim(-0.05, 0.25)
        # ax1.set_ylim(-0.9, 0.4)
        # ax1.set_xlabel("y")
        # ax1.set_ylabel("v")

        # # Prepare line for the solution
        # solution_line, = ax1.plot([], [], color='black', lw=2)
        # critical_manifold, = ax1.plot(og_y_vals, v_vals, color='grey')

        # # Initialize function for FuncAnimation
        # def init():
        #     solution_line.set_data([], [])
        #     critical_manifold.set_data(og_y_vals, v_vals)
            
        #     return solution_line, critical_manifold

        # # Update function for each frame
        # def update(frame):
        #     idx = frame_indices[frame]
        #     solution_line.set_data(actual_sol.y[1][:idx], actual_sol.y[0][:idx])

        #     if actual_sol.t[idx] > 2000:
        #         critical_manifold.set_data(y_vals_2d, v_vals)
            
        #     return solution_line, critical_manifold

        # # Create the animation
        # anim = animation.FuncAnimation(fig, update, frames=len(frame_indices), init_func=init, blit=True, interval=10)

        # # writer = animation.PillowWriter(fps=50)
        # # anim.save("one-spike-phase.gif", writer=writer)
        # plt.show()

        # v_values = np.linspace(-0.9, lower_fold_v_val, points)
        # z_values = np.linspace(0, 0.3, points)
        # # now we have z values, and y values we want the v values (I think)
        # for z_val in z_values: 
        #     for v_val in v_values:
        #         if math.isclose(plot_phi(upper_fold_v_val, z_val), plot_phi(v_val, z_val), abs_tol=0.001):
        #             if not math.isclose(v_val, lower_fold_v_val, abs_tol = 0.0001):
        #                 projected_v_val = v_val

        # projected_fold = np.full(points, projected_v_val)

        # lower_manifold = plot_phi(v_values, z_values)

        # mask = (projected_v_val <= actual_sol.y[0]) & (actual_sol.y[0] <= lower_fold_v_val)
        # phi_val = plot_phi(actual_sol.y[0], actual_sol.y[2])

        # tolerance = 1e-2
        # mask &= np.abs(actual_sol.y[1] - phi_val) < tolerance
        # mask_diff = np.diff(mask.astype(int))

        # transition_indices = np.where(mask_diff != 0)[0]

        # transition_indices = transition_indices - 1
        # transition_indices = transition_indices[transition_indices >= 0]

        # bottom_branch_v = actual_sol.y[0].copy()
        # bottom_branch_z = actual_sol.y[2].copy()

        # bottom_branch_v[~mask] = np.nan
        # bottom_branch_z[~mask] = np.nan

        # bottom_branch_v[transition_indices] = np.nan
        # bottom_branch_z[transition_indices] = np.nan
            
        # if FSN_values:
        #     if stable_inv_manifold:
        #         plt.plot((stable_inv_manifold[0]).y[1], (stable_inv_manifold[0]).y[0], color='purple', label= "Canard solution")
        #     for eq in equilibria:
        #         plt.scatter(plot_zinf(eq), eq, marker='o', color='red', label="True Equilibrium", zorder=10)
        #     plt.plot([0, 0.3], [lower_fold_v_val, lower_fold_v_val], color='green', ls='--', label=r"$F^-$")
        #     plt.plot(z_values, projected_fold, color='black', ls='--', label=r"Projected $F^+$")
        #     plt.plot(bottom_branch_z, bottom_branch_v, color='black')
        #     plt.ylabel('v')
        #     plt.xlabel('z')
        #     plt.ylim(-0.9, -0.3)
        #     plt.xlim(0.0, 0.3)
        #     plt.title(r"Flow on lower stable branch $S^-$ of $M$")
        #     plt.legend(loc="upper right")
        #     plt.grid()
        #     plt.show()