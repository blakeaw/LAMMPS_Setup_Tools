
##########################
# Program to generate tabulated potential/force files for 
# use in LAMMMPS
#
# Adjust relevant parameters and run to generate a tabulated file usable in LAMMPS
# Comes with directory "Potentials" which should be stored in the same directory as 
# this program (e.g. store both users bin directory or standard path directory).  
# Extending:
#   Several Potentials are written and included. You can easily extend by adding new potentials
#   using existing ones as templates. Potential classes are stored in the "Potentials" directory
#   Class file names should be "PotentialClassName.py" 
#   For a potential/force that is a combination of
#   muliple potentials a class can be constructed like the one included
#   "LJ126_COULOMB" 
#	
# Mixing:
#   If your potential has parameters that are the result of
#   mixing others you can either
#   - Precompute the mixed values and give to SetParams
#   - Give the full mixing formula with values to SetParams
#     e.g. Assume LJ126 potential with epsilon1, epsilon2, sigma1, and sigma2.
#          Appy Berthelot mixing by using:
#          potential.SetParams(sqrt(epsilon1*epsilon2), (sigma1+sigma2)*0.5)
#
#  
# Written by Blake Wilson
# email: blake.wilson@utdallas.edu
#
# Please e-mail me if you have comments, suggestions, or experience problems running this.
#
#########################

# import some standard libraries
from math import *


# import potential functions of classes
# Add or adjust imported potential classes as needed
# Potentials directory includes "__init__.py", so is treated
# as module. Must use format :
# from Potentials.PotentClassName import * 
# to import the class of the potential function
# pre-written and included ones are here 
from Potentials.LJ126 import *
from Potentials.LJ96 import *
from Potentials.LJ10050 import *
from Potentials.COULOMB import *
from Potentials.YUKAWA import *
# Example of multipotential class
from Potentials.LJ126_COULOMB import *
from Potentials.LJ10050_YUKAWA import *

# Main

### Parameter setting

#set the bounds of r
# Upper bound
rhigh = 10.0
# Lower Bound
rlow = 1.0
#Number of points to take
npoints = 1000
#Initializes Potential object -- set to desired
potential = LJ126()
#set parameters of potential- adjust values and number of args as needed
potential.SetParams(2.0, 2.0)
#set the name of output file
outname = "something.tab"
# tabulated file type/style line - Adjust as needed 
line = "N " + str(npoints)+" R "+str(rlow)+" "+str(rhigh)
## End of parameter setting

# Open up the file to write to
file_output = open(outname, "w")
# Add some stuff to file to format for LAMMPS tab file
file_output.write("START")
file_output.write("\n")
file_output.write(line)
file_output.write("\n" *2)
# delta r
dr = (rhigh-rlow)/float(npoints-1)

for i in xrange(npoints):
	
#compute current r value
	r = rlow + dr*float(i)
	# Total potential
	U = potential.U(r)
	# Total Force
	F = potential.F(r)
	#write to file
	line = str(i+1)+" "+str(r)+" "+str(U)+" "+str(F)
	file_output.write(line)
	file_output.write("\n")

file_output.close()

#end



