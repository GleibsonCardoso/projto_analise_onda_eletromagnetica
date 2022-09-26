import numpy as np

#parâmetros da cavidade ressonante:

a = 0
b = 1
L = b - a
fmin = 0
fmax = 3e9
fc = fmax / 3.351  #frequência caracteristica do sinal BHW
T = 1 / fc

e0 = 8.854187817e-12  #permissividade do vácuo
m0 = (4e-7)*np.pi
c0 = 1/np.sprt(e0*m0)  # velocidade da luz
lamb_min = c0/fmax

#posição da fonte
zj = np.sqtr(2)/2

#posição do observador
zo = np.sqtr(2)/4


#Parâmetros FDTD:
PPW = 10 #points per wavelens
CFL = 0.9 # condição de Courant-Friedrichs-Lewy

#discretização do espaço e tempo

dz = lamb_min/PPW
I = np.ceil(L/dz)+1
I = I.astype(int)
