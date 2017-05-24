import math

import numpy as np
import matplotlib.pyplot as plt

def analyzer(c, s, cHWP, sHWP, T0):
    # modulator rotates polarization
    c2 =  cHWP*c + sHWP*s
    s2 = -sHWP*c + cHWP*s
    x = c2*c2 - s2*s2 # cos(2(x-a))
    y = 2*c2*s2 # sin(2(x-1))

    # Malus law
    r0 = np.random.random()
    if x > r0*2.0-1.0:
        j = 0 # +1 event
    else:
        j = 1 # -1 event

    # delay time: T0 * sin(2(theta1-x))**4
    l = y*y*y*y*T0*np.random.random()

    return j, l

nsteps = 32
nsamples = 200000
HWP2 = 0
T0 = 1000
W = 1

cHWP2 = math.cos(HWP2/180*math.pi)
sHWP2 = math.sin(HWP2/180*math.pi)
count = np.zeros((2, 2, 2, nsteps))

for ipsi0 in range(nsteps):
    # loop over all different settings of electro-optic modulator
    cHWP1 = math.cos(ipsi0*2.0*math.pi/nsteps)
    sHWP1 = math.sin(ipsi0*2.0*math.pi/nsteps)

    for i in range(nsamples):
        # random polarization for one photon
        angle = np.random.random() * 2.0*math.pi
        c1 = math.cos(angle)
        s1 = math.sin(angle)
        # the other photon is orthogonal
        c2 = -s1
        s2 = c1

        j1, t1 = analyzer(c1, s1, cHWP1, sHWP1, T0)
        j2, t2 = analyzer(c2, s2, cHWP2, sHWP2, T0)

        count[j1, j2, 0, ipsi0] += 1
        if abs(t1-t2) < W:
            count[j1, j2, 1, ipsi0] += 1

# data analysis
tot = np.zeros((2, nsteps))
E1 = np.zeros((2, nsteps))
E2 = np.zeros((2, nsteps))
E12 = np.zeros((2, nsteps))
S = np.zeros((2, nsteps))
r0 = np.zeros(nsteps)
phi = np.zeros(nsteps)

for j in range(nsteps):
    phi[j] = j * 360.0 / nsteps
    r0[j] = -math.cos(2.0*j*2.0*math.pi/nsteps)

    for i in range(2): # temporally correlated or not
        tot[i,j] = np.sum(count[:,:,i,j])
        E12[i,j] = count[0,0,i,j] + count[1,1,i,j] - count[1,0,i,j] - count[0,1,i,j]
        E1[i,j] = count[0,0,i,j] + count[0,1,i,j] - count[1,1,i,j] - count[1,0,i,j]
        E2[i,j] = count[0,0,i,j] + count[1,0,i,j] - count[1,1,i,j] - count[0,1,i,j]

        if tot[i,j] > 0:
            E12[i,j] = E12[i,j] / tot[i,j]
            E1[i,j] = E1[i,j] / tot[i,j]
            E2[i,j] = E2[i,j] / tot[i,j]

        S[i,j] = 3.0 * E12[i,j] - E12[i, j*3 % nsteps]

# plot
#plt.rc('text', usetex=True)
f, axarr = plt.subplots(1, 2, figsize=(15,10))

axarr[0].plot(phi, E12[0,:], label="no time coincidence", linestyle="--", marker="o", color="r")
axarr[0].plot(phi, E12[1,:], label="time coincidence", linestyle="-", marker="o", color="k")
axarr[0].axhline(0, linestyle="--", color="k")
axarr[0].set_xlabel(r"$\varphi$ (degrees)")
axarr[0].set_ylabel(r"$E_{12}(a,b)$")
axarr[0].set_ylim(-1, 1)
axarr[0].legend()

axarr[1].plot(phi, E1[0,:]*E2[0,:], label="no time coincidence", linestyle="--", marker="o", color="r")
axarr[1].plot(phi, E1[1,:]*E2[1,:], label="time coincidence", linestyle="-", marker="o", color="k")
axarr[1].set_xlabel(r"$\varphi$ (degrees)")
axarr[1].set_ylabel(r"$E_1(a,b) \cdot E_2(a,b)$")
axarr[1].set_ylim(-1, 1)
axarr[1].legend()
plt.savefig("eprb_plot.pdf")
plt.show()