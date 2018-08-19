import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from matplotlib.widgets import Slider
from sklearn.preprocessing import normalize
j = 1.0j
#==============================================================================
# preparing the data
#==============================================================================
data = np.loadtxt("002_10C.txt")

freq = data[:,2]
reZ = data[:,0]
imZ = data[:,1]
index = np.arange(0,15)
reZ = np.delete(reZ, index)
imZ = np.delete(imZ, index)
freq = np.delete(freq, index)

omega = 2*np.pi*freq
Z = reZ + j*imZ

#==============================================================================
# defining the function of total impedance
#==============================================================================
def Ztotal(params):
    R1 = params[0]
    R2 = params[1]
    C1 = params[2]
    C2 = params[3]
    Z1 = 1/(j*omega*C1)
    Z2 = R1
    Z3 = 1/(j*omega*C2) + R2
    return (1/Z1 + 1/Z2 + 1/Z3)**-1

#==============================================================================
# defining the function we want to minimize
#==============================================================================
def residual(params):
    logarithmic1abs = np.log(np.sum((np.abs(Z-Ztotal(params))**2)))
    return logarithmic1abs
#==============================================================================
# create the initial guesses
#==============================================================================
params0 = [2.62e7,8.6e5,2.62e-8,2.21e-7]
params00 = [1,1,1,1]
paramparampam = [1.66178211e+06, 5.60598802e+05, 4.24441182e-08, 2.74004036e-07]
initial = [1.0e7, 1.0e6, 5.0e-7, 5.0e-6]
#==============================================================================
# adding the constraints
#==============================================================================
# bounds
# R1bound = (1.0e6,1.0e8)
# R2bound = (1.0e5,1.0e6)
# C1bound = (5.0e-8,1.0e-7)
# C2bound = (5.0e-8,1.0e-6)
# bnds = (R1bound, R2bound, C1bound, C2bound)


# In[89]:


sol = optimize.minimize(residual, params0 , method = "Nelder-Mead" )
print(sol)


# the output is given below

# In[90]:


#print(sol)
def plotting2(params, imZ, reZ):
    labels = ["R1", "R2", "C1", "C2"]
    printing = open("002_10Csol.txt","w")
    for paramparamparampampampampaaaam in np.arange(0,4):
        printout = "{} = {:.2E} \n".format(labels[paramparamparampampampampaaaam],params[paramparamparampampampampaaaam])
        printing.write(printout)
        print(printout)
        #printing.write(labels[paramparamparampampampampaaaam] + " " + str(sol.x[paramparamparampampampampaaaam]) + "\n")
    #==============================================================================
    # PLOTTING
    #==============================================================================
    plt.figure()
    plt.plot(np.real(Ztotal(params)),-1*np.imag(Ztotal(params)), reZ, -imZ, "bo")
    plt.xlabel("Re(Z)")
    plt.ylabel("Im(Z)")
    plt.savefig("002_10C_Nyquist.png", format = "png", dpi = 1200)
    #plt.plot(reZ, -imZ)

    plt.figure()
    plt.plot(np.log(omega), np.log(np.sqrt(reZ**2 + imZ**2)),"bo", np.log(omega), np.log(np.abs(Ztotal(params))), "r-")
    plt.xlabel("log(2*pi*freq)")
    plt.ylabel("log(|Z|)")
    plt.savefig("002_10C_logabsZ.png", format = "png", dpi = 1200)

    plt.figure()
    plt.plot(np.log(omega), np.log(reZ),"bo", np.log(omega), np.log(np.real(Ztotal(params))), "r-")
    plt.xlabel("log(2*pi*freq)")
    plt.ylabel("log(real(|Z|))")
    plt.savefig("002_10C_logrealZ.png", format = "png", dpi = 1200)

    plt.figure()
    plt.plot(np.log(omega), np.log(imZ),"bo", np.log(omega), np.log(-np.imag(Ztotal(params))), "r-")
    plt.xlabel("log(2*pi*freq)")
    plt.ylabel("log(imag|Z|)")
    plt.savefig("002_10C_logimagZ.png", format = "png", dpi = 1200)

    plt.figure()
    plt.plot(np.log(omega), np.arctan(imZ/reZ),"bo", np.log(omega), np.arctan(np.imag(Ztotal(params)/np.real(Ztotal(params)))))
    plt.xlabel("log(2*pi*freq)")
    plt.ylabel("theta")
    plt.savefig("002_10C_theta.png", format = "png", dpi = 1200)


# In[91]:


plotting2(sol.x, imZ, reZ)


# In[78]:


#for something in reZ:
#    print(something)


# In[ ]:
