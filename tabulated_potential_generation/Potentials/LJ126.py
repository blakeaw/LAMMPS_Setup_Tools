# Lennard-Jones 12-6 potential
# has paramters epsilon(energy) and sigma(distance)
class LJ126:

	def __init__(self):
		#Define Params and set some default values
		self.epsilon = 1.0
		self.sigma = 1.0

	def SetParams(self, eps, sig):
		self.epsilon=eps
		self.sigma=sig
		return
	
	def U(self, r):
		return 4.0*self.epsilon*( (self.sigma/r)**12 - (self.sigma/r)**6 )

	def F(self, r):
		return -4.0*self.epsilon*( -12.0*(self.sigma**12/r**13) + 6.0*(self.sigma**6/r**7) )

	
