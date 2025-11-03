#===============
# AUTO Demo cusp
#===============

import sys
sys.path.insert(0, '/Users/mikaeladestradis/Desktop/uni/project/autofiles/auto-07p/python')
import auto # type: ignore
import plotting_reduced as plot

# Load the files (.f90, c.) AUTO command interpreter.
reduced_problem = auto.load('reduced')

try:
    plot.plot_one_param("b.run_type_1", ["Type 1 De-Singularised Bifurcation Diagram", "Type 1 De-Singularised Bifurcation Diagram in Variable z"])
    plot.plot_one_param("b.run_type_2", ["Type 2 De-Singularised Bifurcation Diagram of Equilibria", "Type 2 De-Singularised Bifurcation Diagram of Equilibria in Variable z", "Region at Which We Can Expect Transient Spikes"], zoom_type2=True)
    plot.plot_one_param("b.run_type_3", ["Type 3 De-Singularised Bifurcation Diagram of Equilibria", "Type 3 De-Singularised Bifurcation Diagram of Equilibria in Variable z", "Region at Which We Can Expect Multiple Transient Spikes"], zoom_type3=True)
except:
    run_type_3 = auto.run(reduced_problem)

    auto.save(run_type_3,'run_type_3')

    run_type_2 = auto.run(reduced_problem, c='reduced.2')

    auto.save(run_type_2,'run_type_2')

    #type 1 : to show no transient behaviour is possible!
    run_type_1 = auto.run(reduced_problem, c='reduced.3')

    auto.save(run_type_1,'run_type_1')

    plot.plot_one_param("b.run_type_1", ["Type 1 De-Singularised Bifurcation Diagram", "Type 1 De-Singularised Bifurcation Diagram in Variable z"])
    plot.plot_one_param("b.run_type_2", ["Type 2 De-Singularised Bifurcation Diagram of Equilibria", "Type 2 De-Singularised Bifurcation Diagram of Equilibria in Variable z", "Region at Which We Can Expect Transient Spikes"], zoom=True)
    plot.plot_one_param("b.run_type_3", ["Type 3 De-Singularised Bifurcation Diagram of Equilibria", "Type 3 De-Singularised Bifurcation Diagram of Equilibria in Variable z", "Region at Which We Can Expect Multiple Transient Spikes"], zoom_type3=True)

#clean the directory
auto.wait()
auto.cl()
