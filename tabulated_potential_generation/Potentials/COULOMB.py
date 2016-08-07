# Coulumb potential for electrostatics
# has paramters prefactor(energy*distance)
class COULOMB:

	def __init__(self):
		#Define Params and set some default values
		self.prefactor=1.0

	def SetParams(self, pref):
		self.prefactor=pref
		return

	def U(self, r):
		return self.prefactor/r

	def F(self, r):
		return -self.prefactor/r**2

	
