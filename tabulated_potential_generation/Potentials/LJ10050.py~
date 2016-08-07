# Lennard-Jones 100-50 potential
# has paramters epsilon(energy) and sigma(distance)
class LJ10050:

	def __init__(self):
		#Define Params and set some default values
		self.epsilon = 1.0
		self.sigma = 1.0

	def SetParams(self, eps, sig):
		self.epsilon=eps
		self.sigma=sig
		return

	def U(self, r):
		return 4.0*self.epsilon*( (self.sigma/r)**100 - (self.sigma/r)**50 )

	def F(self, r):
		return -4.0*self.epsilon*( -100.0*(self.sigma**100/r**101) + 50.0*(self.sigma**50/r**51) )

	
