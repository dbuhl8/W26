import matplotlib.pyplot as plt
import numpy as np

fn_grid = ['grid_rfw_plm_minmod_hll.dat', 'grid_rfw_plm_minmod_roe.dat']
fn_slug = ['slug_rfw_plm_minmod_hll_tot.dat', 'slug_rfw_plm_minmod_roe_tot.dat']
sbplt_title = ['PLM + Minmod + Hll', 'PLM + Minmod + ROE']

fig, ax = plt.subplots()

for i in range(1):
    file = open(fn_grid[i],'r')
    nt, nx, xstart, xstop = [float(x) for x in file.readline().split()]
    nt = int(nt)
    nx = int(nx)
    tstop = 0.15

    t = np.zeros(nt)
    x = np.linspace(xstart, xstop, nx+1, True) + (xstop-xstart)/(2.*nx)
    x = x[:nx]

    dens = np.zeros((nt,nx))
    vel = np.zeros_like(dens)
    pres = np.zeros_like(dens)
    gamc = np.zeros_like(dens)
    game = np.zeros_like(dens)
    eint = np.zeros_like(dens)

    file = open(fn_slug[i],'r')

    tidx = 0
    xidx = 0

    for line in file:
        tval, xval, d, v, p, c, e, ei = [float(x) for x in line.split()]
        if not tval in t: 
            if tstop in t:
                break
            tidx += 1
            t[tidx] = tval
        xidx = np.argmax(x>=xval)
        dens[tidx,xidx] = d
        vel[tidx,xidx] = v
        pres[tidx,xidx] = p
        gamc[tidx,xidx] = c
        game[tidx,xidx] = e
        eint[tidx,xidx] = ei

    print('Finished reading file: ',fn_slug[i])

    tidx = np.argmax(t == tstop) - 1
    #print(t)
    #print(tidx)
    ax.plot(x, dens[tidx,:], 'ro-', x, vel[tidx,:], 'go-', x, pres[tidx,:], 'bo-')
    ax.set_title(sbplt_title[i])

plt.savefig('part2_sol.png',dpi=800)

    
