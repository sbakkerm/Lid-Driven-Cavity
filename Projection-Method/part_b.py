import numpy as np
import matplotlib.pyplot as plt

nx = 40 # Points in x
ny = nx # Points in y
xmin = 0 # Minimum coordinate in x
xmax = 1 # Maximum coordinate in x
ymin = 0 # Minimum coordinate in y
ymax = 1 # Maximum coordinate in y
dt = 0.005 # Time step
nstep = 4000 # Number of time steps
Re = 400 # Reynolds Number 
mu = 1/Re # Diffusive coeff
un = 1 # North boundary for U
us = 0 # South boundary for U
ve = 0 # East boundary for V
vw = 0 # West boundary for V
maxit = 1000 # Max iterations for GS/SOR pressure solver
omega = 0.9 # SOR coeff
max_err = 10**(-4) # Convergence criteria for pressure solver

dx = (xmax - xmin)/nx # Grid space in x
dy = (ymax - ymin)/ny # Grid space in y
xs = np.linspace(xmin, xmax, nx+1)
ys = np.linspace(ymin, ymax, ny+1)
x, y = np.meshgrid(xs, ys)


u = np.zeros([nx+1, ny+2]) # Velocity U
v = np.zeros([nx+2, ny+1]) # Velocity V
p = np.ones([nx+2, ny+2]) # Pressure P

pp = np.zeros([nx+1, ny+1]) # Final P
ut = np.zeros([nx+1, ny+2]) # Temp U
vt = np.zeros([nx+2, ny+1]) # Temp V
uu = np.zeros([nx+1, ny+1]) # Final U
vv = np.zeros([nx+1, ny+1]) # Final V

def enforce_boundaries(u, v):
    
    # Enforcing boundaries by interpolation
    u[0:nx+1, 0] = 2*us - u[0:nx+1, 1]
    u[0:nx+1, ny+1] = 2*un - u[0:nx+1, ny]
    v[0, 0:ny+1] = 2*vw - v[1, 0:ny+1]
    v[nx+1, 0:ny+1] = 2*ve - v[nx, 0:ny+1]  
    
    return u, v


def non_vectorized_temp_u(u, v, ut):
    # Temporary u-velocity
    
    for i in range(1, nx):
        for j in range(1, ny+1):
            
            advective = -0.25*(
                ((u[i+1,j] + u[i,j])**2 - (u[i,j] + u[i-1, j])**2)/dx +
                ((u[i,j+1] + u[i,j])*(v[i+1,j] + v[i,j]) - (u[i,j] +u[i,j-1])*(v[i+1,j-1] + v[i,j-1]))/dy)
            
            diffusive = mu*((u[i+1,j] - 2*u[i,j] + u[i-1,j])/dx**2 +
                    (u[i,j+1] - 2*u[i,j] + u[i,j-1])/dy**2)
            
            ut[i,j] = u[i,j] + dt*( advective + diffusive )
            
    return ut

def non_vectorized_temp_v(u, v, vt):
    # Temporary v-velocity

    for i in range(1, nx+1):
        for j in range(1, ny):
            
            advective = -0.25*(
                ((u[i,j+1] + u[i,j])*(v[i+1,j] + v[i,j]) - (u[i-1,j+1] + u[i-1,j])*(v[i,j] + v[i-1,j]))/dx +
                ((v[i,j+1] + v[i,j])**2 - (v[i,j] + v[i,j-1])**2)/dy)
            
            diffusive = mu*((v[i+1,j] -2*v[i,j] + v[i-1,j])/dx**2 +
                    (v[i,j+1] - 2*v[i,j] + v[i,j-1])/dy**2)
            
            vt[i,j] = v[i,j] + dt* ( advective + diffusive )
    
    return vt

def temp_u(u, v, ut):
    # Temporary u-velocity (vectorized)

    advective_u = -0.25*(
                ((u[2:nx+1,1:ny+1] + u[1:nx,1:ny+1])**2 - (u[1:nx,1:ny+1] + u[0:nx-1, 1:ny+1])**2)/dx +
                ((u[1:nx,2:ny+2] + u[1:nx,1:ny+1])*(v[2:nx+1,1:ny+1] + v[1:nx,1:ny+1]) - (u[1:nx,1:ny+1] +u[1:nx,0:ny])*(v[2:nx+1,0:ny] + v[1:nx,0:ny]))/dy)
    
    diffusive_u = mu*((u[2:nx+1,1:ny+1] - 2*u[1:nx,1:ny+1] + u[0:nx-1,1:ny+1])/dx**2 +
                    (u[1:nx,2:ny+2] - 2*u[1:nx,1:ny+1] + u[1:nx,0:ny])/dy**2)
    
    ut[1:nx, 1:ny+1] = u[1:nx, 1:ny+1] + dt*(advective_u + diffusive_u)
    
    return ut
    
def temp_v(u, v, vt):
    # Temporary v-velocity (vectorized)

    advective_v = -0.25*(
                ((u[1:nx+1,2:ny+1] + u[1:nx+1,1:ny])*(v[2:nx+2,1:ny] + v[1:nx+1,1:ny]) - (u[0:nx,2:ny+1] + u[0:nx,1:ny])*(v[1:nx+1,1:ny] + v[0:nx,1:ny]))/dx +
                ((v[1:nx+1,2:ny+1] + v[1:nx+1,1:ny])**2 - (v[1:nx+1,1:ny] + v[1:nx+1,0:ny-1])**2)/dy)
    
    diffusive_v = mu*((v[2:nx+2,1:ny] -2*v[1:nx+1,1:ny] + v[0:nx,1:ny])/dx**2 +
                    (v[1:nx+1,2:ny+1] - 2*v[1:nx+1,1:ny] + v[1:nx+1,0:ny-1])/dy**2)
    
    vt[1:nx+1,1:ny] = v[1:nx+1,1:ny] + dt* ( advective_v + diffusive_v ) 
    
    return vt

def pressure_solver(p, ut, vt):
    
    pt = p.copy()
    
    for it in range(maxit):
        
        pt[0,:] = pt[1,:]
        pt[nx+1,:] = pt[nx,:]
        pt[:,0] = pt[:,1]
        pt[:,ny+1] = pt[:,ny]
        
        for i in range(1, nx+1):
            for j in range(1, ny+1):
                
                # Gauss Seidel
                gs = (0.5/(dx**2 + dy**2)) * ( (dy**2)*(pt[i+1,j] + pt[i-1,j]) +
                                                   (dx**2)*(pt[i,j+1] + pt[i,j-1]) -
                                                   (dx*dy/dt)*(dy*(ut[i,j] - ut[i-1,j]) +
                                                               dx*(vt[i,j] - vt[i,j-1])))
                # SOR (keep as GS with omega = 0)
                pt[i,j] = (1+omega)*gs - omega*pt[i,j]
                
        error = np.abs(pt - p)
        err = np.amax(error)

        p = pt.copy()

        if err < max_err:
            break

    #print(f"n = {n}, error = {err}, iterations = {it}")
    
    return p

"""
MAIN CODE
"""
for n in range(nstep):
    
    u, v = enforce_boundaries(u, v)
    
    # Temporary velocities
    ut = temp_u(u, v, ut)
    vt = temp_v(u, v, vt)
     
    # Solve for Pressure
    p = pressure_solver(p, ut, vt)

    # Correct the velocity
    u[1:nx, 1:ny+1] = ut[1:nx, 1:ny+1] - (dt/dx)*(p[2:nx+1, 1:ny+1] - p[1:nx, 1:ny+1])
    v[1:nx+1, 1:ny] = vt[1:nx+1, 1:ny] - (dt/dy)*(p[1:nx+1, 2:ny+1] - p[1:nx+1, 1:ny])
    
    print(n)
    

# Enforcing boundaries by interpolation
u, v = enforce_boundaries(u, v) 
        
# Restrict the u, v, p grids
uu[0:nx+1, 0:ny+1] = 0.5*(u[0:nx+1, 1:ny+2] + u[0:nx+1, 0:ny+1])
vv[0:nx+1, 0:ny+1] = 0.5*(v[1:nx+2, 0:ny+1] + v[0:nx+1, 0:ny+1])
pp = 0.5*(p[0:nx+1, 0:ny+1] + p[1:nx+2, 1:ny+2])
uu = np.rot90(uu)
vv = np.rot90(vv)
pp = np.rot90(pp)

# Plot final u
plt.figure(dpi=800)
plt.imshow(uu, cmap='turbo', extent = [xmin, xmax, ymin, ymax])
plt.colorbar()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$x$-velcocity: $u$")
plt.show()

# Plot final v
plt.figure(dpi=800)
plt.imshow(vv, cmap='turbo', extent = [xmin, xmax, ymin, ymax])
plt.colorbar()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"$y$-velcocity: $v$")
plt.show()

# Plot final p
plt.figure(dpi=800)
plt.imshow(pp, cmap='coolwarm', extent = [xmin, xmax, ymin, ymax])
plt.colorbar()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"Pressure")
plt.show()              

# Plot streamlines
plt.figure(figsize = (5,5), dpi = 800)
plt.streamplot(x, y, np.flipud(uu), np.flipud(vv), color='black') 
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title(r"Vector Flow Streamlines")
plt.show()       

N = nx + 1


half = int(np.floor(N/2))
vertical = uu[:,half]
horizontal = vv[half,:]

if N%2 == 0:
    
    vertical = (uu[:,half] + uu[:,half-1])/2
    horizontal = (vv[half,:] + vv[half-1,:])/2
    
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
            
    
