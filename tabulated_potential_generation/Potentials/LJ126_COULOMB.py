# Lennard-Jones 12-6 potential
# has paramters epsilon(energy) and sigma(distance)
#  plus
# Coulumb potential for electrostatics
# has paramters prefactor(energy*distance)
class LJ126_COULOMB:

	def __init__(self):
		#Define Params and set some default values
		self.epsilon = 1.0
		self.sigma = 1.0
		self.prefactor=1.0
		

	def SetParams(self, eps, sig, pref):
		self.epsilon=eps
		self.sigma=sig
		self.prefactor=pref
		return

	def U(self, r):
		ULJ = 4.0*self.epsilon*( (self.sigma/r)**12 - (self.sigma/r)**6 )
		UCOUL = self.prefactor/r
		return ULJ+UCOUL

	def F(self, r):
		FLJ = -4.0*self.epsilon*( -12.0*(self.sigma**12/r**13) + 6.0*(self.sigma**6/r**7) )
		FCOUL = -self.prefactor/r**2
		return FLJ+FCOUL

	
