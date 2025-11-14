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

# # Type 1
# try:
#     plot_one_param("b.type_1", "Type 1 Bifurcation")
#     plot_one_param_zoom("b.type_1", "Type 1 Bifurcation in Reasonable Stimulus Range", [-0.05, 0.04])
#     movement_of_equilibria("b.type_1", "Equilibria in Type 1 Parameter Setup", [-0.02, 0.04])
#     plot_frequency("b.type_1", "Frequency Stimulus Plot Under Type 1 Conditions")
# except:
#     # Run and store the result in the Python variable type_1
#     type_1 = auto.run(twodim)

#     auto.save(type_1,'type_1')

#     # save hopf point for two parameter continuation & periodic solutions
#     hb1 = auto.load(type_1('HB1'))

#     # #compute periodic solutions.
#     type_1 = type_1 + auto.run(hb1,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=1000000,DS=-1e-05,DSMAX=1e-03, NTST=50, ITWIST=1)

#     auto.save(type_1,'type_1')

#     plot_one_param("b.type_1", "Type 1 Bifurcation")
#     plot_one_param_zoom("b.type_1", "Type 1 Bifurcation in Reasonable Stimulus Range", [-0.05, 0.04])

#     movement_of_equilibria("b.type_1", "Equilibria in Type 1 Parameter Setup", [-0.02, 0.04])

# # Type 3
# try:
#     plot_one_param("b.type_3", "Type 3 Bifurcation")
# except:
#     #first start with running type 3 (the 'easy' one), here beta_w = -21
#     type_3 = auto.run(twodim, c='twodimnodimens_3')

#     hb3 = auto.load(type_3('HB1'))

#     type_3 = type_3 + auto.run(hb3,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=100000,DS=-1e-05,DSMAX=1e-03)

#     auto.save(type_3,'type_3')

#     plot_one_param("b.type_3", "Type 3 Bifurcation")

# # type 2
# try:
#     plot_one_param("b.type_2", "Type 2 Bifurcation")
#     plot_one_param_zoom("b.type_2", "Type 2 Bifurcation Zoomed in on Hopf point", [0.020, 0.025, -0.5, -0.3], True)
#     movement_of_equilibria("b.type_2", "Equilibria in Type 2 Parameter Setup", [-0.02, 0.04])
#     plot_frequency("b.type_2", "Frequency Stimulus Plot Under Type 2 Conditions")
 
# except:
#     #Now we run with beta w at -13, to show type 2, with the hopf branch as well!
#     type_2 = auto.run(twodim, c='twodimnodimens_2')

#     #save hopf point for two parameter continuation & periodic solutions
#     hb2 = auto.load(type_2('HB1'))

#     print("computing periodic solutions for Type II")
#     #compute periodic solutions.
#     type_2 = type_2 + auto.run(hb2,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=800000,DS=-1e-05,DSMAX=1e-03, NTST=50)

#     auto.save(type_2,'type_2')

#     plot_one_param("b.type_2", "Type 2 Bifurcation")
#     plot_one_param_zoom("b.type_2", "Type 2 Bifurcation Zoomed in on Hopf point", [0.020, 0.025, -0.5, -0.3], True)

#     # last of the one parameter plots:
#     movement_of_equilibria("b.type_2", "Equilibria in Type 2 Parameter Setup", [-0.02, 0.04])

# # two parameter plots     
# try:
#     plot_two_param('b.two_par_cusp', "Two Param Saddle Node Curves")
#     plot_two_param("b.two_par_hopf", "Two Param Hopf Curves")
#     plot_two_param("b.homo", "Homo")
#     plot_two_param("b.two_par_full", "Full Two Parameter Bifurcation Within Reasonable Range")
#     plot_two_param("b.two_par_full", "Two Parameter Zoom on Codimension Two Points", zoom=True, zoom_args=[-10, -8, 0.015, 0.02])
# except:
#     type_1_twop = auto.run(twodim)

#     auto.save(type_1_twop,'type_1_twop')

#     # Set the new start label to the first LP. 
#     lp1 = auto.load(type_1_twop('LP1'))

#     # Continue from this label in two parameters
#     two_par_cusp = auto.run(lp1, ICP=['I', 'beta', 13], NCOL=6, ISW=2, DS = 1e-3, DSMIN= 1e-07, DSMAX= 0.01, EPSL= 1e-09, EPSU = 1e-09, EPSS =1e-06, NTST=100, MXBF=10, IAD=3, IID=3, NMX = 10000, NPR = 10000, UZSTOP = {1: [0.0, 0.04], 2:[-2.1, 0.0]})
    
#     auto.save(two_par_cusp,'two_par_cusp')

#     # Plot the cusp, with Bogdanov Takens
#     plot_two_param('b.two_par_cusp', "Two Param Saddle Node Curves")

#     type_2_twop = auto.run(twodim, c='twodimnodimens_2')

#     hb2 = auto.load(type_2_twop('HB1'))

#     auto.save(type_2_twop,'type_2_twop')

#     two_par_hopf = auto.run(hb2, ICP=[1, 2, 13], NCOL=6, ISW=2, DS = 1e-3, DSMIN= 1e-07, DSMAX= 0.01, EPSL= 1e-09, EPSU = 1e-09, EPSS =1e-06, NTST=100, MXBF=10, IAD=3, IID=3, NMX = 10000000, NPR = 10000000, UZSTOP = {1: [0.0, 1], 2:[-2.5, 2.5]})
#     two_par_hopf = two_par_hopf + auto.run(hb2,ICP=[1, 2, 13], NCOL=6, ISW=2, DS = -1e-3, DSMIN= 1e-07, DSMAX= 0.01, EPSL= 1e-09, EPSU = 1e-09, EPSS =1e-06, NTST=100, MXBF=10, IAD=3, IID=3, NMX = 10000000, NPR = 10000000, UZSTOP = {1: [0.0, 1], 2:[-2.5, 2.5]})

#     auto.save(two_par_hopf,'two_par_hopf')

#     # Plot the hopf curve
#     plot_two_param("b.two_par_hopf", "Two Param Hopf Curves")

#     two_par = two_par_hopf + two_par_cusp

#     bt = auto.load(two_par('BT1'))

#     homo = auto.run(bt, ICP=[1,2,13,11], NUNSTAB=1, NSTAB=2, IEQUIB=1, ITWIST=1, ISTART=1, IREV=[], IFIXED=[13], IPSI=[9,10])

#     auto.save(homo, "homo")

#     plot_two_param("b.homo", "Homo")

#     two_par_full = two_par + homo

#     auto.save(two_par_full, "two_par_full")

#     plot_two_param("b.two_par_full", "Two Par full!")
#     plot_two_param("b.two_par_full", "Two Parameter Zoom on Codimension Two Points", zoom=True, zoom_args=[-10, -8, 0.015, 0.02])

# # plot above the BT
# try: 
#     plot_one_param("b.before_bt", "Bifurcation Plot for Beta Value Less Than BT Point")
#     plot_frequency("b.before_bt", "Frequency Stimulus Plot for Beta Value Less Than BT Point")
# except: 
#     # Run and store the result in the Python variable type_1
#     before_bt = auto.run(twodim, c='before_bt')

#     auto.save(before_bt,'before_bt')

#     # save hopf point for two parameter continuation & periodic solutions
#     hb_before_bt = auto.load(before_bt('HB1'))

#     # #compute periodic solutions.
#     before_bt = before_bt + auto.run(hb_before_bt,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=1000000,DS=-1e-05,DSMAX=1e-03)

#     auto.save(before_bt,'before_bt')

#     plot_one_param("b.before_bt", "Bifurcation Plot for Beta Value Less Than BT Point")
#     plot_frequency("b.before_bt", "Frequency Stimulus Plot for Beta Value Less Than BT Point")

# # plot below the BT
# try: 
#     plot_one_param("b.after_bt", "Bifurcation Plot for Beta Value Larger Than BT Point")
#     plot_frequency("b.after_bt", "Frequency Stimulus Plot for Beta Value Larger Than BT Point")
# except: 
#     after_bt = auto.run(twodim, c='after_bt')

#     auto.save(after_bt,'after_bt')

#     # save hopf point for two parameter continuation & periodic solutions
#     hb1_after_bt = auto.load(after_bt('HB1'))
#     hb2_after_bt = auto.load(after_bt('HB2'))

#     # #compute periodic solutions.
#     after_bt = after_bt + auto.run(hb1_after_bt,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=1000000,DS=-1e-05,DSMAX=1e-03)
#     after_bt = after_bt + auto.run(hb2_after_bt,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=1000000,DS=-1e-05,DSMAX=1e-03)

#     auto.save(after_bt,'after_bt')

#     plot_one_param("b.after_bt", "Bifurcation Plot for Beta Value Larger Than BT Point")
#     plot_frequency("b.after_bt", "Frequency Stimulus Plot for Beta Value Larger Than BT Point")

# # plot below the BT and above the Cusp
# try: 
#     plot_one_param("b.before_cusp", "Bifurcation Plot for Beta Value Less Than Cusp Point")
#     plot_frequency("b.before_cusp", "Frequency Stimulus Plot for Beta Value Less Than Cusp Point")
# except: 
#     before_cusp = auto.run(twodim, c='before_cusp')

#     auto.save(before_cusp,'before_cusp')

#     # save hopf point for two parameter continuation & periodic solutions
#     hb_before_cusp1 = auto.load(before_cusp('HB1'))
#     hb_before_cusp2 = auto.load(before_cusp('HB2'))

#     # #compute periodic solutions.
#     before_cusp = before_cusp + auto.run(hb_before_cusp1,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=1000000,DS=-1e-05,DSMAX=1e-03)
#     before_cusp = before_cusp + auto.run(hb_before_cusp2,IPS=2,ICP=[1,11,12,13],ILP=1,NMX=1000000,DS=-1e-05,DSMAX=1e-03)

#     auto.save(before_cusp,'before_cusp')

#     plot_one_param("b.before_cusp", "Bifurcation Plot for Beta Value Less Than Cusp Point")
#     plot_frequency("b.before_cusp", "Frequency Stimulus Plot for Beta Value Less Than Cusp Point")

# # plot below the cusp
# try: 
#     # either side of the BT
#     plot_one_param("b.after_cusp", "Bifurcation Plot for Beta Value Larger Than Cusp Point")
#     plot_frequency("b.after_cusp", "Frequency Stimulus Plot for Beta Value Larger Than Cusp Point")
# except: 
after_cusp = auto.run(twodim, c='after_cusp')

auto.save(after_cusp,'after_cusp')

# save hopf point for two parameter continuation & periodic solutions
hb_after_cusp = auto.load(after_cusp('HB1'))

# #compute periodic solutions.
# after_cusp = after_cusp + auto.run(hb_after_cusp,IPS=2,ICP=[1,11,12,13],ILP=0,NMX=10000000,DS=1e-04,DSMAX=1e-03)
after_cusp = after_cusp + auto.run(hb_after_cusp,IPS=2,ICP=[1,11,12,13],ILP=0,NMX=10000000,DS=-1e-04,DSMAX=1e-03)

auto.save(after_cusp,'after_cusp')

plot_one_param("b.after_cusp", "Bifurcation Plot for Beta Value Larger Than Cusp Point")
plot_frequency("b.after_cusp", "Frequency Stimulus Plot for Beta Value Larger Than Cusp Point")

# try: 
#     plot_one_param("b.nicer", "More Obvious Type 2")
#     plot_frequency("b.nicer", "Frequency Stimulus Plot for Type 2")
# except:
#     nicer = auto.run(twodim, c='furthertype2')

#     auto.save(nicer,'nicer')

#     # save hopf point for two parameter continuation & periodic solutions
#     nicerhb = auto.load(nicer('HB1'))

#     # #compute periodic solutions.
#     nicer = nicer + auto.run(nicerhb,IPS=2,ICP=[1,11,12,13],ILP=0,NMX=1000000,DS=1e-04,DSMAX=1e-03, NTST=50)

#     auto.save(nicer,'nicer')

#     plot_one_param("b.nicer", "Hopefully works!")

#clean the directory
auto.cl()
auto.wait()
