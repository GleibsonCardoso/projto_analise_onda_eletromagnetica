import numpy as np

#Função de atualização espacial dos campos

#fonte (source)

def s(fc,t):
    a0 = 0.35322222
    a1 = -0.488
    a2 = 0,145
    a3 = -0.01022222
    
    T = 1/fc
    
    soluc = a0 + a1*np.cos(2*np.pi*t/T) + a2*np.cos(2*2*np.pi*t/T) + a3* np.cos(2*3*np.pi*t/T)
    
    return soluc