import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import math
import util
from scipy.integrate import solve_ivp
from systems.system import System
import two_dimensional_system_pars

points = 10000
kv=100
eq_guesses = [-80/kv, -20/kv, -10/kv, 20/kv, 0/kv]
fold_guesses = [10/kv, -40/kv]

w, v = sp.symbols('w v')
I, eps, Vna, gs, Vs, gl, Vl = sp.symbols('I epsilon V_na s V_s l V_l')
gamma_m, m, beta, gamma_b = sp.symbols('gamma_m, m, beta, gamma_b')

winf = (1/2)*(1 + sp.tanh((kv*v - beta)/gamma_b))
dvdt = (I - ((1/2)*(1+sp.tanh((kv*v-m)/gamma_m))*(kv*v - Vna)) - gs*w*(kv*v - Vs) - gl*(kv*v - Vl))/kv
dwdt = eps*(winf-w)*(sp.cosh((kv*v - beta)/(2*gamma_b)))
phi = ((I - (1/2)*(1+sp.tanh((kv*v-m)/gamma_m))*(kv*v - Vna) - gl*(kv*v - Vl)))/(gs*(kv*v - Vs))

params = two_dimensional_system_pars.type_3_start

I_stim = params.get("I")

sympy_params = {sp.symbols(k): v for k, v in params.items()}
exclude = ['I', 'beta', 'epsilon']
sympy_params_exc = {k: v for k, v in params.items() if k not in exclude}

phi_subs = phi.subs(sympy_params)
plot_phi = sp.lambdify(v, phi_subs)
phi_subs_comp = phi.subs(sympy_params_exc)
plot_phi_type = sp.lambdify((v, I, beta), phi_subs_comp)

winf_subs = winf.subs(sympy_params)
plot_winf = sp.lambdify(v, winf_subs)
winf_subs_comp = winf.subs(sympy_params_exc)
plot_winf_type = sp.lambdify((v, I, beta), winf_subs_comp)

#equilibria occur when: g(V, phi(V)) = 0 (ie when the nullclines cross)
eq_curve = dvdt.subs(w, winf)
eq_curve = eq_curve.subs(sympy_params)
equilibria = util.find_equilibria(eq_curve, v, eq_guesses, I_stim)

# finding the folds & plotting full nullcline solutions
fold_points = util.find_fold_points(dvdt, phi, w, v, sympy_params, fold_guesses, I_stim)   

# set up to plot
v_vals = np.linspace(-99/kv, 50/kv, points)
phi_vals = plot_phi(v_vals)
g_vals = plot_winf(v_vals)

# set up the plot of flow when epsilon is turned on
f_func = sp.lambdify((v, w, I), dvdt.subs(sympy_params_exc))
g_func = sp.lambdify((v, w, beta, eps), dwdt.subs(sympy_params_exc))

def twodim(t, X, I, b, eps):
    v, w = X

    return[float(f_func(v, w, I)), float(g_func(v, w, b, eps))]


def setup_nullcline_plot(title, alt_title):
    plt.title(title)
    plt.xlabel("w")
    plt.ylabel("v")
    plt.xlim(-0.1, 0.3)
    plt.ylim(-0.9, 0.5)
    plt.legend(loc='upper right')
    try:
        plt.savefig("pics/" + title)
    except: 
        plt.savefig("pics/" + alt_title)
    plt.close()

def shoot_up(init, stimulus, beta):
    # compute the shoot up: straight up to the manifold, so the w value is unchanging. 
    start_v = init[0]
    start_w = init[1]
    shoot_up = np.linspace(start_v, 0.5, points)
    for v in shoot_up:
        if math.isclose(plot_phi_type(v, stimulus, beta), start_w, abs_tol=0.0005):
            return (v)
    return []

def reduced_flow_upper(init, stimulus, beta):
    start_v = init[0]
    folds = util.find_fold_points(dvdt, phi, w, v, sympy_params_exc, fold_guesses, stimulus)
    v_values = np.linspace(start_v, folds[0], points)
    w_values = plot_phi_type(v_values, stimulus, beta)

    return [v_values, w_values]

def jump_point(init, stimulus, beta):
    start_v = init[0]
    start_w = init[1]
    drop_off = np.linspace(-0.9, start_v, points)
    for v in drop_off:
        if math.isclose(plot_phi_type(v, stimulus, beta), start_w, abs_tol=0.001):
            return (v)
    return []

def reduced_flow_lower(init, stimulus, beta):
    start_v = init[0]
    equilibria_curve = dvdt.subs(w, winf)
    equilibria_curve = equilibria_curve.subs(sympy_params_exc)
    equilibria_curve = equilibria_curve.subs('beta', beta)
    eq = util.find_equilibria(equilibria_curve, v, eq_guesses, stimulus)
    
    folds = util.find_fold_points(dvdt, phi, w, v, sympy_params_exc, fold_guesses, stimulus)

    if eq[0] < folds[1]:
        v_values = np.linspace(start_v, eq[0], points)
        w_values = plot_phi_type(v_values, stimulus, beta)
    else:
        v_values = np.linspace(start_v, folds[1], points)
        w_values = plot_phi_type(v_values, stimulus, beta)

    return [v_values, w_values]

def shoot_back(init, stimulus, beta):
    start_v = init[0]
    start_w = init[1]
    shoot_up = np.linspace(start_v+0.1, 0.5, points)
    for v in shoot_up:
        if math.isclose(plot_phi_type(v, stimulus, beta), start_w, abs_tol=0.001):
            return (v)
    return []

def plot_time_traces(solutions, title):
    colours=['red', 'blue']
    labels = ["epsilon = 0.015", "epsilon = 0.0015"]
    for i in range(len(solutions)):
        sol=solutions[i]
        plt.plot(sol.t, sol.y[0], color=colours[i], label=labels[i])
    plt.xlabel("time (T)")
    plt.xlim(0, 2000)
    plt.ylabel("v")
    plt.title(title)
    plt.legend(loc="lower right")
    plt.savefig("pics/" + title)
    plt.close()

def compute_solution(initial_conditions, I, b, eps):
    t_span = [0, 2000]
    points_dt= int(6000/(0.001))
    t_eval = np.linspace(t_span[0], t_span[1], points_dt)

    sol = solve_ivp(
        twodim,
        t_span,
        initial_conditions,
        method="LSODA",
        t_eval=t_eval,
        args=(I, b, eps)
    )

    return sol

# V = kv*v, w = w
class NonDimensionalTwoDim(System):
    def run():
        plt.plot(phi_vals, v_vals, color='black', label="Slow Manifold")
        plt.plot(g_vals, v_vals, color='green', label="Nullcline")
        plt.plot(plot_phi(fold_points[0]), fold_points[0], marker='o', color='darkviolet', label=r'$F^{+}$')
        plt.plot(plot_phi(fold_points[1]), fold_points[1], marker='o', color='red', label=r'$F^{-}$')
        if len(equilibria) == 1:
            plt.plot(plot_phi(equilibria[0]), equilibria[0],  marker='o', color='blue', label='equilibrium')
        else: 
            for eq in equilibria:
                plt.plot(plot_phi(eq), eq,  marker='o', color='red', label='equilibrium')
        plt.xlabel('w')
        plt.ylabel('v')
        plt.title("Nullcline for Type 1")
        plt.xlim(-0.1, 0.35)
        plt.ylim(-92/kv,50/kv)
        plt.legend()
        plt.grid()
        plt.savefig("General Nullcine")
        plt.close()

        def mini_nullclines():
            plt.xlabel("w")
            plt.ylabel("v")
            plt.xlim(-0.1, 0.3)
            plt.ylim(-0.9, 0.5)

        #side by side plots of how the nullcline changes - talk to Martin to check this is alright to put in/how I can make sure it's ref if need be.
        plt.plot(plot_phi_type(v_vals, 80/(20), 0), v_vals, color='black')
        plt.plot(plot_winf_type(v_vals, 80/(20), 0), v_vals, color='green')
        mini_nullclines()
        plt.savefig("pics/2 dimensional system/side_by_side/80,0")
        plt.close()

        plt.plot(plot_phi_type(v_vals, 40/(20), 0), v_vals, color='black')
        plt.plot(plot_winf_type(v_vals, 40/(20), 0), v_vals, color='green')
        mini_nullclines()
        plt.savefig("pics/2 dimensional system/side_by_side/40,0")
        plt.close()

        plt.plot(plot_phi_type(v_vals, 0, 0), v_vals, color='black')
        plt.plot(plot_winf_type(v_vals, 0, 0), v_vals, color='green')
        mini_nullclines()
        plt.savefig("pics/2 dimensional system/side_by_side/0,0")
        plt.close()

        plt.plot(plot_phi_type(v_vals, 80/(20), -13), v_vals, color='black')
        plt.plot(plot_winf_type(v_vals, 80/(20), -13), v_vals, color='green')   
        mini_nullclines()    
        plt.savefig("pics/2 dimensional system/side_by_side/80,-13")
        plt.close()
        
        plt.plot(plot_phi_type(v_vals, 40/(20), -13), v_vals, color='black')
        plt.plot(plot_winf_type(v_vals, 40/(20), -13), v_vals, color='green')
        mini_nullclines()
        plt.savefig("pics/2 dimensional system/side_by_side/40,-13")
        plt.close()

        plt.plot(plot_phi_type(v_vals, 0, -13), v_vals, color='black')
        plt.plot(plot_winf_type(v_vals, 0, -13), v_vals, color='green')    
        mini_nullclines()  
        plt.savefig("pics/2 dimensional system/side_by_side/0,-13")
        plt.close()

        # plots for 'type 3'
        type_3_eq_curve = dvdt.subs(w, winf)
        type_3_eq_curve = type_3_eq_curve.subs(sympy_params_exc)
        type_3_eq_curve = type_3_eq_curve.subs(beta, -21)
        type3_kickoff = util.find_equilibria(type_3_eq_curve, v, eq_guesses, I_stim)
        initial_conditions = [type3_kickoff[0], plot_winf_type(type3_kickoff[0], 0, -21)]
        eq_20 = util.find_equilibria(type_3_eq_curve, v, eq_guesses, 20/20)
        eq_50 = util.find_equilibria(type_3_eq_curve, v, eq_guesses, 50/20)
        eq_80 = util.find_equilibria(type_3_eq_curve, v, eq_guesses, 80/20)

        plt.plot(plot_winf_type(v_vals, 0, -21), v_vals, color='blue') 
        plt.plot(plot_phi_type(v_vals, 0, -21), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 20/(20), -21), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 50/(20), -21), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 80/(20), -21), v_vals, color='black')
        plt.scatter(initial_conditions[1], initial_conditions[0], color='red', label='Equilibria', zorder=10, s=18)
        plt.scatter(plot_winf_type(eq_20[0], 20/20, -21), eq_20[0], color='red', zorder=10, s=18)
        plt.scatter(plot_winf_type(eq_50[0], 50/20, -21), eq_50[0], color='red', zorder=10, s=18)
        plt.scatter(plot_winf_type(eq_80[0], 80/20, -21), eq_80[0], color='red', zorder=10, s=18)

        setup_nullcline_plot("Nullcines for Type 3, with Varying Stimulus Value", "Nullcines for Type 3, with Varying Stimulus Value")

        final_20 = shoot_up(initial_conditions, 20/(20), -21)
        final_50 = shoot_up(initial_conditions, 50/(20), -21)
        final_80 = shoot_up(initial_conditions, 80/(20), -21)

        upper_50 = reduced_flow_upper([final_50, initial_conditions[1]], 50/(20), -21)
        upper_80 = reduced_flow_upper([final_80, initial_conditions[1]], 80/(20), -21)

        dropped_50 = jump_point([upper_50[0][-1], upper_50[1][-1]], 50/(20), -21)
        dropped_80 = jump_point([upper_80[0][-1], upper_80[1][-1]], 80/(20), -21)

        lower_50 = reduced_flow_lower([dropped_50, upper_50[1][-1]], 50/(20), -21)
        lower_80 = reduced_flow_lower([dropped_80, upper_80[1][-1]], 80/(20), -21)

        shoot_back_80 = shoot_back([lower_80[0][-1], lower_80[1][-1]], 80/(20), -21)

        solution_20_high_eps = compute_solution(initial_conditions, 20/20, -21, params.get("epsilon"))
        solution_50_high_eps = compute_solution(initial_conditions, 50/20, -21, params.get("epsilon"))
        solution_80_high_eps = compute_solution(initial_conditions, 80/20, -21, params.get("epsilon"))
        solution_20_low_eps = compute_solution(initial_conditions, 20/20, -21, (params.get("epsilon"))/10)
        solution_50_low_eps = compute_solution(initial_conditions, 50/20, -21, (params.get("epsilon"))/10)
        solution_80_low_eps = compute_solution(initial_conditions, 80/20, -21, (params.get("epsilon"))/10)

        # Plot some time traces: 
        plot_time_traces([solution_20_high_eps, solution_20_low_eps], "Time trace for stimulus 20 and beta -21")
        plot_time_traces([solution_50_high_eps, solution_50_low_eps], "Time trace for stimulus 50 and beta -21")
        plot_time_traces([solution_80_high_eps, solution_80_low_eps], "Time trace for stimulus 80 and beta -21")

        # plot of GSPT style (erpsilon limit) flow for stimulus 80$\mu A/cm^2$ (or 0.04 in non-dimensionalised parameters)
        plt.plot(plot_winf_type(v_vals, 0, -21), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='orange', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 80/(20), -21), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_80], color="green")
        plt.plot(upper_80[1], upper_80[0], color='green')
        plt.plot([upper_80[1][-1], upper_80[1][-1]], [upper_80[0][-1], dropped_80], color='green')
        plt.plot(lower_80[1], lower_80[0], color='green', label="Singular Limit Solution")
        plt.plot([lower_80[1][-1], lower_80[1][-1]], [lower_80[0][-1], shoot_back_80], color='green')
        plt.scatter(lower_80[1][-1], lower_80[0][-1], color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_80_high_eps.y[1], solution_80_high_eps.y[0], color='red', label="Solution for epsilon 0.015")
        plt.plot(solution_80_low_eps.y[1], solution_80_low_eps.y[0], color='blue', label="Solution for epsilon 0.0015")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=80$\mu A/cm^2$ and beta=-21", "Nullclines and GSPT flow for I=80 and beta=-21")

        # this is interesting behaviour, let's zoom in a bit. plot the fold with the equilibria and that...
        new_params = sympy_params.copy()
        new_params.update({beta: -21, eps: 0.0015, I: 80/20})
        folds_int = util.find_fold_points(dvdt, phi, w, v, new_params, fold_guesses, 80/(20))   
        eq_int = util.find_equilibria(type_3_eq_curve, v, eq_guesses, 80/(20))

        plt.plot(plot_phi_type(v_vals, 80/(20), -21), v_vals, color='black')
        plt.plot(lower_80[1], lower_80[0], color='green')
        plt.plot([lower_80[1][-1], lower_80[1][-1]], [lower_80[0][-1], shoot_back_80], color='green', label="Singular Limit Solution")
        plt.scatter(plot_winf_type(eq_int[0], 80/20, -21), eq_int[0], color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_80_low_eps.y[1], solution_80_low_eps.y[0], color='purple', label="Solution for epsilon 0.0015")
        plt.scatter(plot_phi_type(folds_int[1], 80/20, -21), folds_int[1], color='orange', label="Fold point", zorder=12, alpha=0.9)
        plt.xlim(0.025, 0.075)
        plt.ylim(-0.6, -0.2)
        plt.xlabel("w")
        plt.ylabel("v")
        plt.legend(loc='upper right')
        plt.title("Interactions Between the Lower Fold and Equilibrium")
        plt.savefig("Interactions Between the Lower Fold and Equilibrium")
        plt.close()


        # plot of GSPT style (erpsilon limit) flow for stimulus 50$\mu A/cm^2$ (or 0.025 in non-dimensionalised parameters)
        plt.plot(plot_winf_type(v_vals, 0, -21), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='orange', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 50/(20), -21), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_50], color="green")
        plt.plot(upper_50[1], upper_50[0], color='green')
        plt.plot([upper_50[1][-1], upper_50[1][-1]], [upper_50[0][-1], dropped_50], color='green')
        plt.plot(lower_50[1], lower_50[0], color='green', label="Singular Limit Solution")
        plt.scatter(lower_50[1][-1], lower_50[0][-1], color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_50_high_eps.y[1], solution_50_high_eps.y[0], color='red', label="Solution for epsilon 0.015")
        plt.plot(solution_50_low_eps.y[1], solution_50_low_eps.y[0], color='blue', label="Solution for epsilon 0.0015")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=50$\mu A/cm^2$ and beta=-21", "Nullclines and GSPT flow for I=50 and beta=-21")

        # plot of GSPT style (erpsilon limit) flow for stimulus 20$\mu A/cm^2$ (or 0.01 in non-dimensionalised parameters)
        plt.plot(plot_winf_type(v_vals, 0, -21), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='orange', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 20/(20), -21), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_20], color="green", label="Singular Limit Solution")
        plt.scatter(initial_conditions[1], final_20, color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_20_high_eps.y[1], solution_20_high_eps.y[0], color='red', label="Solution for epsilon 0.015")
        plt.plot(solution_20_low_eps.y[1], solution_20_low_eps.y[0], color='blue', label="Solution for epsilon 0.0015")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=20$\mu A/cm^2$ and beta=-21", "Nullclines and GSPT flow for I=20 and beta=-21")        

        # Now for type 2: 
        type_2_eq_curve = dvdt.subs(w, winf)
        type_2_eq_curve = type_2_eq_curve.subs(sympy_params_exc)
        type_2_eq_curve = type_2_eq_curve.subs(beta, -13)
        type2_kickoff = util.find_equilibria(type_2_eq_curve, v, eq_guesses, I_stim)
        initial_conditions = [type2_kickoff[0], plot_winf_type(type2_kickoff[0], 0, -13)]
        eq_20 = util.find_equilibria(type_2_eq_curve, v, eq_guesses, 20/20)
        eq_50 = util.find_equilibria(type_2_eq_curve, v, eq_guesses, 50/20)
        eq_80 = util.find_equilibria(type_2_eq_curve, v, eq_guesses, 80/20)

        plt.plot(plot_winf_type(v_vals, 0, -13), v_vals, color='blue') 
        plt.plot(plot_phi_type(v_vals, 0, -13), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 20/(20), -13), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 50/(20), -13), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 80/(20), -13), v_vals, color='black')
        plt.scatter(initial_conditions[1], initial_conditions[0], color='red', label='Equilibria', zorder=10, s=18)
        plt.scatter(plot_winf_type(eq_20[0], 20/20, -13), eq_20[0], color='red', zorder=10, s=18)
        plt.scatter(plot_winf_type(eq_50[0], 50/20, -13), eq_50[0], color='red', zorder=10, s=18)
        plt.scatter(plot_winf_type(eq_80[0], 80/20, -13), eq_80[0], color='red', zorder=10, s=18)
        setup_nullcline_plot("Nullcines for Type 2, with Varying Stimulus Value", "Nullcines for Type 2, with Varying Stimulus Value")

        final_20 = shoot_up(initial_conditions, 20/(20), -13)
        final_50 = shoot_up(initial_conditions, 50/(20), -13)
        final_80 = shoot_up(initial_conditions, 80/(20), -13)

        upper_50 = reduced_flow_upper([final_50, initial_conditions[1]], 50/(20), -13)
        upper_80 = reduced_flow_upper([final_80, initial_conditions[1]], 80/(20), -13)

        dropped_50 = jump_point([upper_50[0][-1], upper_50[1][-1]], 50/(20), -13)
        dropped_80 = jump_point([upper_80[0][-1], upper_80[1][-1]], 80/(20), -13)

        lower_50 = reduced_flow_lower([dropped_50, upper_50[1][-1]], 50/(20), -13)
        lower_80 = reduced_flow_lower([dropped_80, upper_80[1][-1]], 80/(20), -13)

        shoot_back_50 = shoot_back([lower_50[0][-1], lower_50[1][-1]], 50/(20), -13)
        shoot_back_80 = shoot_back([lower_80[0][-1], lower_80[1][-1]], 80/(20), -13)

        solution_20_high_eps = compute_solution(initial_conditions, 20/20, -13, params.get('epsilon'))
        solution_50_high_eps = compute_solution(initial_conditions, 50/20, -13, params.get('epsilon'))
        solution_80_high_eps = compute_solution(initial_conditions, 80/20, -13, params.get('epsilon'))
        solution_20_low_eps = compute_solution(initial_conditions, 20/20, -13, params.get('epsilon')*0.1)
        solution_50_low_eps = compute_solution(initial_conditions, 50/20, -13, params.get('epsilon')*0.1)
        solution_80_low_eps = compute_solution(initial_conditions, 80/20, -13, params.get('epsilon')*0.1)

        # Plot some time traces: 
        plot_time_traces([solution_20_high_eps, solution_20_low_eps], "Time trace for stimulus 20 and beta -13")
        plot_time_traces([solution_50_high_eps, solution_50_low_eps], "Time trace for stimulus 50 and beta -13")
        plot_time_traces([solution_80_high_eps, solution_80_low_eps], "Time trace for stimulus 80 and beta -13")

        plt.plot(plot_winf_type(v_vals, 0, -13), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='orange', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 80/(20), -13), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_80], color="green")
        plt.plot(upper_80[1], upper_80[0], color='green')
        plt.plot([upper_80[1][-1], upper_80[1][-1]], [upper_80[0][-1], dropped_80], color='green')
        plt.plot(lower_80[1], lower_80[0], color='green')
        plt.plot([lower_80[1][-1], lower_80[1][-1]], [lower_80[0][-1], shoot_back_80], color='green', label="Singular Limit Solution")
        plt.scatter(lower_80[1][-1], lower_80[0][-1], color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_80_high_eps.y[1], solution_80_high_eps.y[0], color='red', label="Solution for epsilon 0.015")
        plt.plot(solution_80_low_eps.y[1], solution_80_low_eps.y[0], color='blue', label="Solution for epsilon 0.0015")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=80$\mu A/cm^2$ and beta=-13", "Nullclines and GSPT flow for I=80 and beta=-13")

        plt.plot(plot_winf_type(v_vals, 0, -13), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='orange', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 50/(20), -13), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_50], color="green")
        plt.plot(upper_50[1], upper_50[0], color='green')
        plt.plot([upper_50[1][-1], upper_50[1][-1]], [upper_50[0][-1], dropped_50], color='green')
        plt.plot(lower_50[1], lower_50[0], color='green')
        plt.plot([lower_50[1][-1], lower_50[1][-1]], [lower_50[0][-1], shoot_back_50], color='green', label="Singular Limit Solution")
        plt.scatter(lower_50[1][-1], lower_50[0][-1], color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_50_high_eps.y[1], solution_50_high_eps.y[0], color='red', label="Solution for epsilon 0.015")
        plt.plot(solution_50_low_eps.y[1], solution_50_low_eps.y[0], color='blue', label="Solution for epsilon 0.0015")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=50$\mu A/cm^2$ and beta=-13", "Nullclines and GSPT flow for I=50 and beta=-13")

        plt.plot(plot_winf_type(v_vals, 0, -13), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='orange', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 20/(20), -13), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_20], color="green", label="Singular Limit Solution")
        plt.scatter(initial_conditions[1], final_20, color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_20_high_eps.y[1], solution_20_high_eps.y[0], color='red', label="Solution for epsilon 0.015")
        plt.plot(solution_20_low_eps.y[1], solution_20_low_eps.y[0], color='blue', label="Solution for epsilon 0.0015")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=20$\mu A/cm^2$ and beta=-13", "Nullclines and GSPT flow for I=20 and beta=-13")

        # Finally for type 1: 
        type_1_eq_curve = dvdt.subs(w, winf)
        type_1_eq_curve = type_1_eq_curve.subs(sympy_params_exc)
        type_1_eq_curve = type_1_eq_curve.subs(beta, 0)
        type1_kickoff = util.find_equilibria(type_1_eq_curve, v, eq_guesses, I_stim)
        initial_conditions = [type1_kickoff[0], plot_winf_type(type1_kickoff[0], 0, 0)]
        eq_20 = util.find_equilibria(type_1_eq_curve, v, eq_guesses, 20/20)
        eq_40 = util.find_equilibria(type_1_eq_curve, v, eq_guesses, 40/20)
        eq_80 = util.find_equilibria(type_1_eq_curve, v, eq_guesses, 80/20)

        plt.plot(plot_winf_type(v_vals, 0, 0), v_vals, color='blue') 
        plt.plot(plot_phi_type(v_vals, 0, 0), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 20/(20), 0), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 40/(20), 0), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 80/(20), 0), v_vals, color='black')
        plt.scatter(initial_conditions[1], initial_conditions[0], color='red', label="Equilibria", zorder=10,  s=18)
        plt.scatter(plot_winf_type(type1_kickoff[1], 0, 0), type1_kickoff[1], color='red', zorder=10,  s=18)
        plt.scatter(plot_winf_type(type1_kickoff[2], 0, 0), type1_kickoff[2], color='red', zorder=10,  s=18)
        for eq in eq_20:
            plt.scatter(plot_winf_type(eq, 20/(20), 0), eq, color='red', zorder=10, s=18)
        for eq in eq_40:
            plt.scatter(plot_winf_type(eq, 40/(20), 0), eq, color='red', zorder=10, s=18)
        for eq in eq_80:
            plt.scatter(plot_winf_type(eq, 80/(20), 0), eq, color='red', zorder=10, s=18)
        setup_nullcline_plot("Nullcines for Type 1, with Varying Stimulus Value", "Nullcines for Type 1, with Varying Stimulus Value")

        final_20 = shoot_up(initial_conditions, 20/(20), 0)
        final_40 = shoot_up(initial_conditions, 40/(20), 0)
        final_80 = shoot_up(initial_conditions, 80/(20), 0)

        upper_40 = reduced_flow_upper([final_40, initial_conditions[1]], 40/(20), 0)
        upper_80 = reduced_flow_upper([final_80, initial_conditions[1]], 80/(20), 0)

        dropped_40 = jump_point([upper_40[0][-1], upper_40[1][-1]], 40/(20), 0)
        dropped_80 = jump_point([upper_80[0][-1], upper_80[1][-1]], 80/(20), 0)

        lower_40 = reduced_flow_lower([dropped_40, upper_40[1][-1]], 40/(20), 0)
        lower_80 = reduced_flow_lower([dropped_80, upper_80[1][-1]], 80/(20), 0)

        shoot_back_40 = shoot_back([lower_40[0][-1], lower_40[1][-1]], 40/(20), 0)
        shoot_back_80 = shoot_back([lower_80[0][-1], lower_80[1][-1]], 80/(20), 0)

        solution_20_high_eps = compute_solution(initial_conditions, 20/20, 0, params.get("epsilon"))
        solution_40_high_eps = compute_solution(initial_conditions, 40/20, 0, params.get("epsilon"))
        solution_80_high_eps = compute_solution(initial_conditions, 80/20, 0, params.get("epsilon"))
        solution_20_low_eps = compute_solution(initial_conditions, 20/20, 0, params.get("epsilon")*0.1)
        solution_40_low_eps = compute_solution(initial_conditions, 40/20, 0, params.get("epsilon")*0.1)
        solution_80_low_eps = compute_solution(initial_conditions, 80/20, 0, params.get("epsilon")*0.1)

        # Plot some time traces: 
        plot_time_traces([solution_20_high_eps, solution_20_low_eps], "Time trace for stimulus 20 and beta 0")
        plot_time_traces([solution_40_high_eps, solution_40_low_eps], "Time trace for stimulus 40 and beta 0")
        plot_time_traces([solution_80_high_eps, solution_80_low_eps], "Time trace for stimulus 80 and beta 0")

        plt.plot(plot_winf_type(v_vals, 0, 0), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='orange', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 80/(20), 0), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_80], color="green")
        plt.plot(upper_80[1], upper_80[0], color='green')
        plt.plot([upper_80[1][-1], upper_80[1][-1]], [upper_80[0][-1], dropped_80], color='green')
        plt.plot(lower_80[1], lower_80[0], color='green')
        plt.plot([lower_80[1][-1], lower_80[1][-1]], [lower_80[0][-1], shoot_back_80], color='green', label="Singular Limit Solution")
        plt.plot(solution_80_high_eps.y[1], solution_80_high_eps.y[0], color='red', label="Solution for epsilon 0.015")
        plt.plot(solution_80_low_eps.y[1], solution_80_low_eps.y[0], color='blue', label="Solution for epsilon 0.0015")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=80$\mu A/cm^2$ and beta=0", "Nullclines and GSPT flow for I=80 and beta=0")

        plt.plot(plot_winf_type(v_vals, 0, 0), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='orange', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 40/(20), 0), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_40], color="green")
        plt.plot(upper_40[1], upper_40[0], color='green')
        plt.plot([upper_40[1][-1], upper_40[1][-1]], [upper_40[0][-1], dropped_40], color='green')
        plt.plot(lower_40[1], lower_40[0], color='green')
        plt.plot([lower_40[1][-1], lower_40[1][-1]], [lower_40[0][-1], shoot_back_40], color='green', label="Singular Limit Solution")
        plt.plot(solution_40_high_eps.y[1], solution_40_high_eps.y[0], color='red', label="Solution for epsilon 0.015")
        plt.plot(solution_40_low_eps.y[1], solution_40_low_eps.y[0], color='blue', label="Solution for epsilon 0.0015")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=40$\mu A/cm^2$ and beta=0", "Nullclines and GSPT flow for I=40 and beta=0")

        plt.plot(plot_winf_type(v_vals, 0, 0), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='orange', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 20/(20), 0), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_20], color="green", label="Singular Limit Solution")
        plt.scatter(initial_conditions[1], final_20, color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_20_high_eps.y[1], solution_20_high_eps.y[0], color='red', label="Solution for epsilon 0.015")
        plt.plot(solution_20_low_eps.y[1], solution_20_low_eps.y[0], color='blue', label="Solution for epsilon 0.0015")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=20 $\mu A/cm^2$ and beta=0", "Nullclines and GSPT flow for I=20 and beta=0")

        # Now I'd also like some stuff near beta=-10!
        type_1_eq_curve = dvdt.subs(w, winf)
        type_1_eq_curve = type_1_eq_curve.subs(sympy_params_exc)
        type_1_eq_curve = type_1_eq_curve.subs(beta, -8)
        type1_kickoff = util.find_equilibria(type_1_eq_curve, v, eq_guesses, I_stim)
        initial_conditions = [type1_kickoff[0], plot_winf_type(type1_kickoff[0], 0, -8)]

        plt.plot(plot_winf_type(v_vals, 0, -8), v_vals, color='blue') 
        plt.plot(plot_phi_type(v_vals, 0,-8), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 20/(20),-8), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 50/(20),-8), v_vals, color='black')
        plt.plot(plot_phi_type(v_vals, 80/(20),-8), v_vals, color='black')
        plt.scatter(initial_conditions[1], initial_conditions[0], color='red', label='Initial Condition', zorder=10)
        setup_nullcline_plot("Nullcines for Type ?, with Varying Stimulus Value", "Nullcines for Type ?, with Varying Stimulus Value")

        final_20 = shoot_up(initial_conditions, 20/(20),-8)
        final_50 = shoot_up(initial_conditions, 50/(20),-8)
        final_80 = shoot_up(initial_conditions, 80/(20),-8)

        upper_50 = reduced_flow_upper([final_50, initial_conditions[1]], 50/(20),-8)
        upper_80 = reduced_flow_upper([final_80, initial_conditions[1]], 80/(20),-8)

        dropped_50 = jump_point([upper_50[0][-1], upper_50[1][-1]], 50/(20),-8)
        dropped_80 = jump_point([upper_80[0][-1], upper_80[1][-1]], 80/(20),-8)

        lower_50 = reduced_flow_lower([dropped_50, upper_50[1][-1]], 50/(20),-8)
        lower_80 = reduced_flow_lower([dropped_80, upper_80[1][-1]], 80/(20),-8)

        shoot_back_50 = shoot_back([lower_50[0][-1], lower_50[1][-1]], 50/(20),-8)
        shoot_back_80 = shoot_back([lower_80[0][-1], lower_80[1][-1]], 80/(20),-8)

        solution_20_high_eps = compute_solution(initial_conditions, 20/20,-8, params.get('epsilon'))
        solution_50_high_eps = compute_solution(initial_conditions, 50/20,-8, params.get('epsilon'))
        solution_80_high_eps = compute_solution(initial_conditions, 80/20,-8, params.get('epsilon'))

        # Plot some time traces: 
        plot_time_traces([solution_20_high_eps], "Time trace for stimulus 20 and beta -8")
        plot_time_traces([solution_50_high_eps], "Time trace for stimulus 50 and beta -8")
        plot_time_traces([solution_80_high_eps], "Time trace for stimulus 80 and beta -8")

        plt.plot(plot_winf_type(v_vals, 0,-8), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='red', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 80/(20),-8), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_80], color="green")
        plt.plot(upper_80[1], upper_80[0], color='green')
        plt.plot([upper_80[1][-1], upper_80[1][-1]], [upper_80[0][-1], dropped_80], color='green')
        plt.plot(lower_80[1], lower_80[0], color='green')
        plt.plot([lower_80[1][-1], lower_80[1][-1]], [lower_80[0][-1], shoot_back_80], color='green', label="Singular Limit Solution")
        plt.scatter(lower_80[1][-1], lower_80[0][-1], color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_80_high_eps.y[1], solution_80_high_eps.y[0], color='red', label="Solution to full problem")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=80$\mu A/cm^2$ and beta=0", "Nullclines and GSPT flow for I=80 and beta=-8")

        plt.plot(plot_winf_type(v_vals, 0,-8), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='red', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 50/(20),-8), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_50], color="green", label="Singular Limit Solution")
        plt.plot(upper_50[1], upper_50[0], color='green')
        plt.plot([upper_50[1][-1], upper_50[1][-1]], [upper_50[0][-1], dropped_50], color='green')
        plt.plot(lower_50[1], lower_50[0], color='green')
        plt.plot([lower_50[1][-1], lower_50[1][-1]], [lower_50[0][-1], shoot_back_50], color='green')
        plt.scatter(lower_50[1][-1], lower_50[0][-1], color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_50_high_eps.y[1], solution_50_high_eps.y[0], color='red', label="Solution to full problem")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=50$\mu A/cm^2$ and beta=0", "Nullclines and GSPT flow for I=50 and beta=-8")

        plt.plot(plot_winf_type(v_vals, 0,-8), v_vals, color='lightblue') 
        plt.scatter(initial_conditions[1], initial_conditions[0], color='red', label='Initial Condition', zorder=10)
        plt.plot(plot_phi_type(v_vals, 20/(20),-8), v_vals, color='black')
        plt.plot([initial_conditions[1], initial_conditions[1]], [initial_conditions[0], final_20], color="green", label="Singular Limit Solution")
        plt.scatter(initial_conditions[1], final_20, color='blue', label='Equilibrium', zorder=10)
        plt.plot(solution_20_high_eps.y[1], solution_20_high_eps.y[0], color='red', label="Solution to full problem")
        setup_nullcline_plot(r"Nullclines and GSPT flow for I=20 $\mu A/cm^2$ and beta=0", "Nullclines and GSPT flow for I=20 and beta=-8")

    def move():
        return
