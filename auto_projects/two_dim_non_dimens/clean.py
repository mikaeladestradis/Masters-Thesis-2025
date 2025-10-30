#===============
# AUTO Demo cusp
#===============

import sys
sys.path.insert(0, '/Users/mikaeladestradis/Desktop/uni/project/autofiles/auto-07p/python')
from auto import * # type: ignore

print("\n***Clean the directory***")

# delete("type_1")
# delete("type_2")
# delete("type_3")
delete("two_par")
delete("two_par_hopf")
delete("two_par_cusp")
clean()