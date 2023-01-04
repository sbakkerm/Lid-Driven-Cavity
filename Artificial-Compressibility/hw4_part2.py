import numpy as np
import matplotlib.pyplot as plt

U = 1
Re = 400 #100
H = 1 # domain
N = 51 # 16, 32, 64, 128
h = H /(N-1) # delta x = delta y
dt = 0.001 
coeff = 0.4 
beta = 1/coeff
P = 1

xs = np.linspace(0, H, N)
ys = np.linspace(0, H, N)
x, y = np.meshgrid(xs, ys)
N = xs.shape[0]

u = np.zeros((N, N))
v = np.zeros((N, N))
p = P*np.ones((N, N))

# y's, x's (j, i or rows, cols)
u[-1, :] = U # initial condition

a = 1
b = N-1


def update_uij(u, v, p, i, j):
    
    uij = u[j, i]
    D_u = (-4*uij + u[j, i+1] + u[j, i-1] + u[j+1, i] + u[j-1, i]) / (Re*h**2)
    C_u = (u[j, i+1]**2 - u[j, i-1]**2 + u[j+1, i]*v[j+1, i] - u[j-1, i]*v[j-1, i]) / (2*h)
    P_u = -(p[j, i+1] - p[j, i-1])/(2*h)
    
    uij_new = uij + dt*(P_u + D_u - C_u)
    
    return uij_new

def update_uij_vectorized(u, v, p):
    
    uij = u[a:b, a:b]
    
    D_u = (-4*uij + u[a:b, a+1:b+1] + 
           u[a:b, a-1:b-1] +
           u[a+1:b+1, a:b] +
           u[a-1:b-1, a:b]) / (Re*h**2)
    
    C_u = (u[a:b, a+1:b+1]**2 - u[a:b, a-1:b-1]**2 +
           u[a+1:b+1, a:b]*v[a+1:b+1, a:b] -
           u[a-1:b-1, a:b]*v[a-1:b-1, a:b]) / (2*h)
    
    P_u = -(p[a:b, a+1:b+1] - p[a:b, a-1:b-1]) / (2*h)
    
    u_new = uij + dt*(P_u + D_u - C_u)
    
    return u_new
    
 
def update_vij(u, v, p, i, j):
    
    vij = v[j, i]
    D_v = (-4*vij + v[j, i+1] + v[j, i-1] + v[j+1, i] + v[j-1, i]) / (Re*h**2)
    C_v = (u[j, i+1]*v[j, i+1] - u[j, i-1]*v[j, i-1] + v[j+1, i]**2 - v[j-1, i]**2) / (2*h)
    P_v = -(p[j+1, i] - p[j-1, i])/(2*h)
    
    vij_new = vij + dt*(P_v + D_v - C_v)
    
    return vij_new

def update_vij_vectorized(u, v, p):
    
    vij = v[a:b, a:b]
    
    D_v = (-4*vij + v[a:b, a+1:b+1] +
           v[a:b, a-1:b-1] +
           v[a+1:b+1, a:b] +
           v[a-1:b-1, a:b]) / (Re*h**2)
    
    C_v = (u[a:b, a+1:b+1]*v[a:b, a+1:b+1] - 
           u[a:b, a-1:b-1]*v[a:b, a-1:b-1] +
           v[a+1:b+1, a:b]**2 - v[a-1:b-1, a:b]**2) / (2*h)
    
    P_v = -(p[a+1:b+1, a:b] - p[a-1:b-1, a:b]) / (2*h)
    
    v_new = vij + dt*(P_v + D_v - C_v)
    
    return v_new

def update_pij(u, v, p, i, j):

    pij = p[j, i] 
    
    dudx = (u[j, i+1] - u[j, i-1])/(2*h)
    dvdx = (v[j+1, i] - v[j-1, i])/(2*h)

        
    pij_new = pij - dt*(dudx + dvdx)/beta
    
    return pij_new

def update_pij_vectorized(u, v, p):
    
    pij = p[a:b, a:b]
    
    dudx = (u[a:b, a+1:b+1] - u[a:b, a-1:b-1]) / (2*h)
    
    dvdx = (v[a+1:b+1, a:b] - v[a-1:b-1, a:b]) / (2*h)
    
    p_new = pij - dt*(dudx + dvdx)/beta

    return p_new    
    


def timeStep(u, v, p):
    
    u_new = np.zeros((N, N))
    v_new = np.zeros((N, N))
    p_new = np.zeros((N, N))
    

    u_new[-1, :] = U
    #u_new[0, :] = U
    
    #errors = np.zeros((N,N))
    error = 0
    
    # Comment out below for non-vectorized version
    """
    for j in range(1,N-1):
        
        for i in range(1,N-1):
            
            #u_new[j, i] = update_uij(u, v, p, i, j)
            #v_new[j, i] = update_vij(u, v, p, i, j)
            p_new[j, i] = update_pij(u, v, p, i, j)
    """
    
    u_new[a:b, a:b] = update_uij_vectorized(u, v, p)
    v_new[a:b, a:b] = update_vij_vectorized(u, v, p)
    p_new[a:b, a:b] = update_pij_vectorized(u, v, p)
    
    p_new[:,0] = p_new[:,1]
    p_new[:,N-1] = p_new[:,N-2]
    p_new[0,:] = p_new[1,:]
    p_new[N-1,:] = p_new[N-2,:]
    
    
    errors = (u_new - u + v_new - v)/dt  # dudt + dvdt
    error = np.sqrt(np.sum( np.sum(errors**2, axis=0)))/(N**2)

    return (u_new, v_new, p_new, error/N**2)
    


error = 1
threshold = 1e-10
iteration = 0
while error > threshold:
    
    u_new, v_new, p_new, error = timeStep(u, v, p)

    u, u_new = u_new, u
    v, v_new, = v_new, v
    p, p_new = p_new, p
    
    iteration += 1
    
    if iteration%1000 == 0:
        print(iteration, error)

   
plt.figure(dpi=800)
plt.imshow(u, cmap = 'jet', extent = [0, 1, 0, 1], origin = 'lower')
plt.colorbar()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$x$-velcocity: $u$")
plt.show()  

plt.figure(dpi = 800)
plt.imshow(v, cmap = 'jet', extent = [0, 1, 0, 1], origin = 'lower')
plt.colorbar()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$y$-velcocity: $v$")
plt.show()  

plt.figure(dpi = 800)
plt.imshow(p, cmap = 'jet', extent = [0, 1, 0, 1], origin = 'lower')
plt.colorbar()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"Pressure")
plt.show()         

plt.figure(figsize = (5,5), dpi = 800)
plt.streamplot(x, y, u, v, color='black') 
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"Vector Flow Streamlines")
plt.show()       

half = int(np.floor(N/2))
vertical = u[:,half]
horizontal = v[half,:]

if N%2 == 0:
    
    vertical = (u[:,half] + u[:,half-1])/2
    horizontal = (v[half,:] + v[half-1,:])/2
    
# vertical centerline at num half
#"""
plt.figure(dpi = 800)
plt.plot(ys, vertical, 'r')
plt.plot(ys, vertical, 'k.')
plt.title(r"$u$ along the the vertical centerline $x=0.5$")
plt.xlabel(r"$y$")
plt.ylabel(r"$x$-velocity: $u$")
plt.show()

print(f"u min is {min(vertical)}")

plt.figure(dpi = 800)
plt.plot(xs, horizontal, 'r')
plt.plot(xs, horizontal, 'k.')
plt.title(r"$v$ along the the horizontal centerline $y=0.5$")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$-velocity: $v$")
plt.show()
#"""
print(f"v min is {min(horizontal)}")
print(f"v max is {max(horizontal)}")

title_u = "Re=" + str(Re) + "_N=" + str(N) + "_u.txt"
title_v = "Re=" + str(Re) + "_N=" + str(N) + "_v.txt"
title_p = "Re=" + str(Re) + "_N=" + str(N) + "_p.txt"

np.savetxt(title_u, u)
np.savetxt(title_v, v)
np.savetxt(title_p, p)

































