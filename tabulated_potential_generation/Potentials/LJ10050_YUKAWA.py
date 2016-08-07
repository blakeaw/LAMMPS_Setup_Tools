# Lennard-Jones 100-50 potential plus Yukawa style potential for screened electrostatics
# has paramters LJ: epsilon(energy) and sigma(distance); Yukawa: prefactor(energy) and kappa (1/distance)
# Can be used as a model for charged colloids in electrolytic solution
class LJ10050_YUKAWA:

	def __init__(self):
		#Define Params and set some default values
		self.epsilon = 1.0
		self.sigma = 1.0
		self.prefactor=1.0
		self.kappa = 1.0

	def SetParams(self, eps, sig, pref, kap):
		self.epsilon=eps
		self.sigma=sig
		self.prefactor=pref
		self.kappa = kap
		return

	def U(self, r):
		ULJ= 4.0*self.epsilon*( (self.sigma/r)**100 - (self.sigma/r)**50 )
		UYK= (self.prefactor/self.kappa)*exp(-self.kappa*r)/r
		return ULJ+UYK

	def F(self, r):
		FLJ= -4.0*self.epsilon*( -100.0*(self.sigma**100/r**101) + 50.0*(self.sigma**50/r**51) )
		FYK = (self.prefactor)*exp(-self.kappa*r)/r + (self.prefactor/self.kappa)*exp(-self.kappa*r)/r**2
		return FLJ+FYK
	
