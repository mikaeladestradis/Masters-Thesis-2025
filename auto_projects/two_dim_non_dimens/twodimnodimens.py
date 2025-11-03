#===============
# AUTO Demo cusp
#===============

import sys
sys.path.insert(0, '/Users/mikaeladestradis/Desktop/uni/project/autofiles/auto-07p/python')
import auto # type: ignore
from plotting import *

# Load the files .f90 and c.xxx into the AUTO
# command interpreter.
twodim = auto.load('twodimnodimens')

# Type 1
try:
    plot_one_param("b.type_1", "Type 1 Bifurcation")
    plot_one_param_zoom("b.type_1", "Type 1 Bifurcation in Reasonable Stimulus Range", [-0.05, 0.04])
    movement_of_equilibria("b.type_1", "Equilibria in Type 1 Parameter Setup", [-0.02, 0.04])
except:
    # Run and store the result in the Python variable type_1
    type_1 = auto.run(twodim)

    auto.save(type_1,'type_1')

    # save hopf point for two parameter continuation & periodic solutions
    hb1 = auto.load(type_1('HB1'))

    # #compute periodic solutions.
    type_1 = type_1 + auto.run(hb1,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=1000000,DS=-1e-05,DSMAX=1e-03)

    auto.save(type_1,'type_1')

    plot_one_param("b.type_1", "Type 1 Bifurcation")
    plot_one_param_zoom("b.type_1", "Type 1 Bifurcation in Reasonable Stimulus Range", [-0.05, 0.04])

    movement_of_equilibria("b.type_1", "Equilibria in Type 1 Parameter Setup", [-0.02, 0.04])

# Type 3
try:
    plot_one_param("b.type_3", "Type 3 Bifurcation")
except:
    #first start with running type 3 (the 'easy' one), here beta_w = -21
    type_3 = auto.run(twodim, c='twodimnodimens_3')

    hb3 = auto.load(type_3('HB1'))

    type_3 = type_3 + auto.run(hb3,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=100000,DS=-1e-05,DSMAX=1e-03)

    auto.save(type_3,'type_3')

    plot_one_param("b.type_3", "Type 3 Bifurcation")

try:
    plot_one_param("b.type_2", "Type 2 Bifurcation")
    plot_one_param_zoom("b.type_2", "Type 2 Bifurcation Zoomed in on Hopf point", [0.020, 0.025, -0.5, -0.3], True)
    movement_of_equilibria("b.type_2", "Equilibria in Type 2 Parameter Setup", [-0.02, 0.04])
 
except:
    #Now we run with beta w at -13, to show type 2, with the hopf branch as well!
    type_2 = auto.run(twodim, c='twodimnodimens_2')

    #save hopf point for two parameter continuation & periodic solutions
    hb2 = auto.load(type_2('HB1'))

    print("computing periodic solutions for Type II")
    #compute periodic solutions.
    type_2 = type_2 + auto.run(hb2,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=800000,DS=-1e-05,DSMAX=1e-03, NTST=50)

    auto.save(type_2,'type_2')

    plot_one_param("b.type_2", "Type 2 Bifurcation")
    plot_one_param_zoom("b.type_2", "Type 2 Bifurcation Zoomed in on Hopf point", [0.020, 0.025, -0.5, -0.3], True)

    # last of the one parameter plots:
    movement_of_equilibria("b.type_2", "Equilibria in Type 2 Parameter Setup", [-0.02, 0.04])
        
try:
    plot_two_param('b.two_par_cusp', "Two Param Saddle Node Curves")
    plot_two_param("b.two_par_hopf", "Two Param Hopf Curves")
    plot_two_param("b.two_par", "Full Two Parameter Bifurcation")
    plot_two_param("b.homo", "Homo")
    plot_two_param("b.two_par_full", "Two Par full!")
except:
    type_1_twop = auto.run(twodim)

    auto.save(type_1_twop,'type_1_twop')

    # Set the new start label to the first LP. 
    lp1 = auto.load(type_1_twop('LP1'))

    # Continue from this label in two parameters
    two_par_cusp = auto.run(lp1, ICP=['I', 'beta', 13], NCOL=6, ISW=2, DS = 1e-3, DSMIN= 1e-07, DSMAX= 0.01, EPSL= 1e-09, EPSU = 1e-09, EPSS =1e-06, NTST=100, MXBF=10, IAD=3, IID=3, NMX = 10000, NPR = 10000, UZSTOP = {1: [0.0, 0.04], 2:[-2.1, 0.0]})
    
    auto.save(two_par_cusp,'two_par_cusp')

    # Plot the cusp, with Bogdanov Takens
    plot_two_param('b.two_par_cusp', "Two Param Saddle Node Curves")

    type_2_twop = auto.run(twodim, c='twodimnodimens_2')

    hb2 = auto.load(type_2_twop('HB1'))

    auto.save(type_2_twop,'type_2_twop')

    two_par_hopf = auto.run(hb2, ICP=[1, 2, 13], NCOL=6, ISW=2, DS = 1e-3, DSMIN= 1e-07, DSMAX= 0.01, EPSL= 1e-09, EPSU = 1e-09, EPSS =1e-06, NTST=100, MXBF=10, IAD=3, IID=3, NMX = 10000000, NPR = 10000000, UZSTOP = {1: [0.0, 1], 2:[-2.5, 2.5]})
    two_par_hopf = two_par_hopf + auto.run(hb2,ICP=[1, 2, 13], NCOL=6, ISW=2, DS = -1e-3, DSMIN= 1e-07, DSMAX= 0.01, EPSL= 1e-09, EPSU = 1e-09, EPSS =1e-06, NTST=100, MXBF=10, IAD=3, IID=3, NMX = 10000000, NPR = 10000000, UZSTOP = {1: [0.0, 1], 2:[-2.5, 2.5]})

    auto.save(two_par_hopf,'two_par_hopf')

    # Plot the hopf curve
    plot_two_param("b.two_par_hopf", "Two Param Hopf Curves")

    two_par = two_par_hopf + two_par_cusp

    bt = auto.load(two_par('BT1'))

    homo = auto.run(bt, NUNSTAB=1, NSTAB=2, IEQUIB=1, ITWIST=1, ISTART=1, IREV=[], IFIXED=[13], IPSI=[9,10])

    auto.save(homo, "homo")

    two_par_full = two_par + homo

    auto.save(two_par_full, "two_par_full")

    plot_two_param("b.two_par_full", "Two Par full!")

    at_bt = auto.run(bt, c="twodimnodimens")

    auto.save(at_bt, "at_bt")

    plot_one_param("b.at_bt", "first attempt- at bt point?")
    # todo: change the colour of the homoclinic branch
    # make a zoomed in picture that is closer to the BT point, and have three bifurcation diagrmas
    #  what happens nearby - 1) the homoclinc side, 2) the BT point itself and 3) the hopf side. 

#clean the directory
auto.cl()
auto.wait()
