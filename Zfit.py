import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import pylab
j = 1.0j


#creating a function of the total impedance

#==============================================================================
#preparing the data
#==============================================================================
data = np.loadtxt("006_30C.txt")

freq = data[:,2]
reZ = data[:,0]
imZ = data[:,1]
index = np.arange(0,15)
reZ = np.delete(reZ, index)
imZ = np.delete(imZ, index)
freq = np.delete(freq, index)

omega = 2*np.pi*freq
Z = reZ + j*imZ



def Ztotal(params):
    R1 = params[0]
    R2 = params[1]
    C1 = params[2]
    C2 = params[3]
    Z1 = 1/(j*omega*C1)
    Z2 = R1
    Z3 = 1/(j*omega*C2) + R2
    return (1/Z1 + 1/Z2 + 1/Z3)**-1

def constraint1(params):
    return params[0] - 80000
def constraint2(params):
    return 120000 - params[0]


#def residual(params):
#    return np.log(np.sum((np.abs(Z) - np.abs(Ztotal(params)))**2))

def residual(params):
    #logarithmic errors
    first = np.log(np.sum((np.abs(Z) - np.abs(Ztotal(params)))**2))
    #impedance errors


params0 = [1.62e6,2.86e5,2.62e-8,2.21e-7]
params00 = [1,1,1,1]
paramparampam = [1.66178211e+06, 5.60598802e+05, 4.24441182e-08, 2.74004036e-07]

con1 = {"type":"ineq", "fun" : constraint1}
con2 = {"type":"ineq", "fun" : constraint2}
cons = [con1, con2]
sol = optimize.minimize(residual, paramparampam, method = "Nelder-Mead" )
print(sol)


labels = ["R1", "R2", "C1", "C2"]
printing = open("sol.txt","w")
for paramparamparampampampampaaaam in np.arange(0,4):
    printing.write(labels[paramparamparampampampampaaaam] + " " + str(sol.x[paramparamparampampampampaaaam]) + "\n")
#==============================================================================
# PLOTTING
#==============================================================================
plt.subplot(2,2,1)
plt.plot(np.real(Ztotal(sol.x)),-1*np.imag(Ztotal(sol.x)),"ro", reZ, -imZ)
plt.xlabel("Re(Z)")
plt.ylabel("Im(Z)")
#plt.plot(reZ, -imZ)
plt.subplot(2,2,2)
plt.plot(np.log(omega), np.log(np.sqrt(reZ**2 + imZ**2)),"bo", np.log(omega), np.log(np.abs(Ztotal(sol.x))), "r-")
plt.xlabel("log(2*pi*freq)")
plt.ylabel("log(|Z|)")

plt.subplot(2,2,3)

plt.plot(np.log(omega), np.log(reZ),"bo", np.log(omega), np.log(np.abs(np.real(Ztotal(sol.x)))), "r-")
plt.xlabel("log(2*pi*freq)")
plt.ylabel("log(real(|Z|))")

plt.subplot(2,2,4)

plt.plot(np.log(omega), np.log(reZ),"bo", np.log(omega), np.log(np.abs(np.imag(Ztotal(sol.x)))), "r-")
plt.xlabel("log(2*pi*freq)")
plt.ylabel("log(imag|Z|)")
plt.savefig("plotted.png", format = "png", dpi = 1200)
