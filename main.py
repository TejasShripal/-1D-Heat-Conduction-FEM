import numpy as np
from pandas import DataFrame

#function to build a tridiagonal matrix
def tridiagonal_matrix(lower, middle, upper):
    n = len(middle)
    matrix = np.zeros((n, n))

    # Fill the main diagonal
    matrix[np.arange(n), np.arange(n)] = middle

    # Fill the lower diagonal
    matrix[np.arange(1, n), np.arange(n-1)] = lower

    # Fill the upper diagonal
    matrix[np.arange(n-1), np.arange(1, n)] = upper

    return matrix

def r(i):
    return dr * (i + 1)

#MODIFY DATA HERE ---------------------------------------------------------------------------------------------/
l = 1
n = 10                                           #Number of elements
t0 = 293                                         #Initial temperature       
tinf = 293                                       #Surrounding temperature
K = 45                                           #Conductivity constant   
del_t = 0.01                                     #Timestep
t_max = 1                                        #Timeframe
rhow = 7750                                      #Density 
cp = 420                                         #Specific Heat Capacity
diff = K / (rhow * cp)                           #Diffusivity constant
B = diff * del_t                                 #Diffusivity constant * timestamp constant
dr = 1/n                                         #Space Step 
c = diff * del_t / (dr**2)                       #Fourier Number
h = 25.32                                        #Convective heat coefficient
q = 63694                                        #Source heat generation
t = 0.00                                         #initial time variable
#--------------------------------------------------------------------------------------------------------------/

#Building Mass Matrix-----------------------------------------------------------------------------------------/
#lower M1
lower = np.zeros(n-1)
for i in range(1,len(lower)+1): # 1 -> 9 since r(i) is from 2 -> 10, refer r(i) function above
    if i == len(lower):
        lower[-1] = dr * (r(i-1) + r(i))/12
    else:
        lower[i-1] = r(i) * dr / 3
#print(lower)

#upper M1
upper = np.zeros(n-1)
for i in range(0,len(upper)): # 
    upper[i] = dr * (r(i+1) + r(i))/12
#print(upper)

#middle M1

middle = np.zeros(n)
for i in range(len(middle)):
    if i == 0:
        middle[i] = ((3 * r(i)) + r(i+1))*dr/12
    elif i == 9:
         middle[i] = ((3 * r(i)) + r(i-1))*dr/12
    else:
        middle[i] = 2 * r(i) * dr / 3
#print(middle)

M1 = tridiagonal_matrix(lower, middle, upper)
print(DataFrame(M1))
#-----------------------------------------------------------------------------------------/

#Finding the inverse of the Mass Matrix
Minv = np.linalg.inv(M1)
print(DataFrame(Minv))



#M2 (Stiffness Matrix + Mass Matrix)

#lower
lower = np.zeros(n-1)
for i in range(1,len(lower)+1): # 1 -> 9 since r(i) is from 2 -> 10, refer r(i) function above
    if i == len(lower):
        lower[-1] = (dr * (r(i-1) + r(i))/12) - (B * (2 * r(i-1) + r(i)) / 6)
    else:
        lower[i-1] = (-1 * (r(i) + 2 * r(i-1)) * B / 6) + (r(i) * dr / 3)

#upper
upper = np.zeros(n-1)
for i in range(0,len(upper)): # 
    if i == 0:
        upper[0] = (B * (r(i) + 2 * r(i+1)) / 6) + (r(i) + r(i+1) * dr / 12)
    else:
        upper[i] = (-1 * B * ( 2 * r(i+1) - r(i)) / 6) + (r(i) + r(i+1) * dr / 12)

#middle
middle = np.zeros(n)
for i in range(0, len(middle)):
    for i in range(len(middle)):
        if i == 0:
            middle[i] = (B* ((2*r(i)) + r(i+1)) / 6) + ((3 * r(i)) + r(i+1))*dr/12
        elif i == 9:
            middle[i] = ((-1 * B * ((2*r(i)) + r(i-1)) / 6) + (B*h)) + ((3 * r(i)) + r(i-1))*dr/12
        else:
            middle[i] = (B * dr / 3) + (2*r(i) * dr / 3)

M2 = tridiagonal_matrix(lower, middle, upper)
print("M22-----------------------")
print(DataFrame(M2))

#-----------------------------------------------------------------------------------------/

#Source Matrix

SM = (B/K) * M1
print(DataFrame(SM))

#Boundary Condition Vector
BC = np.zeros(n)
BC[-1] = (h * tinf  * del_t) #+ (h * tinf * del_t / K)

#Current Time step temperature vector

Tk = t0 * np.ones(n)

#Next time step temperature vector
Tk1 = np.zeros(n)

#Source Vector
S = q * np.ones(n)

#saving results to path
file_path = f'results/{0}.txt'
np.savetxt(file_path, Tk)

#Initializing Loop variables
i = 1
t = 0.01
#Loop to solve for X in a system of linear equations in the form AX = B
while(t < t_max):
    if i == 99:
        break
    C = np.dot(M2,Tk) + np.dot(SM,S)+ BC
    #T[i+1] = np.dot(Minv, C)
    Tk = np.linalg.solve(M1,C)
    file_path = f'results/{i}.txt'
    np.savetxt(file_path, Tk)
    i += 1
    t += del_t

#-----------------------------------------------------------------------------------------/
    

















    


