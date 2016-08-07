# Lennard-Jones 9-6 potential
# has paramters epsilon(energy) and sigma(distance)
class LJ96:

	def __init__(self):
		#Define Params and set some default values
		self.epsilon = 1.0
		self.sigma = 1.0

	def SetParams(self, eps, sig):
		self.epsilon=eps
		self.sigma=sig
		return

	def U(self, r):
		return 4.0*self.epsilon*( (self.sigma/r)**9 - (self.sigma/r)**6 )

	def F(self, r):
		return -4.0*self.epsilon*( -9.0*(self.sigma**9/r**10) + 6.0*(self.sigma**6/r**7) )

	
