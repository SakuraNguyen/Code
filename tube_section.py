import numpy as np

r = 1.
theta_min = 45.
theta_max = 180. - theta_min
npoints = 19 # INTEGER
zl = 1.
nl = 4 # INTEGER

arc_increm = (theta_max - theta_min)/(npoints-1)

out = np.zeros((npoints*nl,3))
for j in range(0,nl):
    for i in range(0,npoints):
        print (j,i)
        pt_x = r*np.cos((theta_min + i*arc_increm)/180. * np.pi)
        pt_y = r*np.sin((theta_min + i*arc_increm)/180. * np.pi)
        pt_z = j*zl 
        out[i+j*npoints] = np.array([pt_x,pt_y,pt_z])
    
print (out)
