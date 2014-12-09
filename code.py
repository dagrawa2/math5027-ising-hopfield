from __future__ import division
from pylab import *

# if you import randint from pylab,
# then randint(a,b) does not include b
# but if you import randint from random,
# then randint(a,b) does include b

class ising():
	def __init__(self,N,T,Train):
		self.N = N  # size of 1D states
		self.L = N-1  # number of time slices per period
		self.T = T  # temperature
		for i in range(len(Train)):
			for j in range(N):
				Train[i][j] = 2*Train[i][j]-1  # change 0s to -1s
		W = zeros((N,N))  # Hopfield weights
		for i in range(N):
			for j in range(N):
				for k in range(len(Train)):
					W[i,j] += Train[k][i]*Train[k][j]
		W = W/len(Train)
		self.J = zeros((N,N))  # spatial couplings
		self.h = zeros((N,N))  # time couplings
		for i in range(N):
			for j in range(N):
				self.J[i,j] = W[(i*j+j)%N,(i*j+j+i+1)%N]
				self.h[i,j] = W[(i*j+j)%N,(i*j+2*j)%N]
		self.state = zeros((N,N-1))  # overall 2D state (one period)

	def initial(self,ini):  # create 2D cylinder
		for i in range(self.N):
			ini[i] = 2*ini[i]-1
		self.state[:,0] = ini  # fix initial 1D state
		for i in range(self.N):
			for j in range(1,self.L):
				self.state[i,j] = 2*randint(0,2)-1  # randomize

	def run(self, num=1):
		for step in range(num):
			z = randint(self.N, self.N*self.L)  # select site
			(i,j) = (z%self.N,z//self.N)
			s = self.state[i,j]
			WF = self.J[j,(i-1)%self.N]*self.state[(i-1)%self.N,j] \
+ self.J[j,i%self.N]*self.state[(i+1)%self.N,j] \
+ self.h[j-1,i]*self.state[i,j-1]
			DeltaE = WF*2*s  # change in energy
			if DeltaE < 0:
				self.state[i,j] = -s
			elif exp(-DeltaE/self.T) > rand():
				self.state[i,j] = -s

	def get(self,steps,freq):
		num = 0
		S = [0]*self.N
		rounds = steps//freq
		for i in range(rounds):
			self.run(freq)
			num += 1
			for j in range(self.L):  # unpermute states
				for i in range(self.N):
					S[(i*j+i)%self.N] += self.state[i,j]
		for i in range(self.N):
			S[i] = S[i]/(num*self.L)  # average over time slices
		return S


train = [ \
[1] + [0]*10, \
[0]*4 + [1,1] + [0]*3 + [1,1]]  # training set
initial = [0]*4 + [1]*7

I = ising(11,0.1,train)

S = initial
for i in range(10):  # number of periods
	I.initial(S)
	I.run(20000)
	S = I.get(20000,10)
	S = [int(sign(S[i])+1)//2 for i in range(len(S))]

print S
