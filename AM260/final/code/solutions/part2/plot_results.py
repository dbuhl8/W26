import matplotlib.pyplot as plt
import numpy as np

pdir = '../../slug_scripts/part2_'
gpre = '/grid_rfw_plm_minmod_'
ppre = '/slug_rfw_plm_minmod_'
ppath = ['hll', 'roe']
ppsuf = '_tot'
psuf = '.dat'
fn_grid = [pdir+x+gpre+x+psuf for x in ppath]
fn_slug = [pdir+x+ppre+x+ppsuf+psuf for x in ppath]
debugging = True
if (debugging):
    fn_grid.append('../../slug_scripts/part2_fog/grid_rfw_fog_minmod_hll.dat')
    fn_slug.append('../../slug_scripts/part2_fog/slug_rfw_fog_minmod_hll_tot.dat')
sbplt_title = ['PLM + Minmod + Hll', 'PLM + Minmod + ROE', 'FOG + HLL Baseline']
#fn_grid = ['grid_rfw_plm_minmod_hll.dat', 'grid_rfw_plm_minmod_roe.dat']
#fn_slug = ['slug_rfw_plm_minmod_hll_tot.dat', 'slug_rfw_plm_minmod_roe_tot.dat']
#sbplt_title = ['PLM + Minmod + Hll', 'PLM + Minmod + ROE']
nr = 3
fig, ax = plt.subplots(nrows=nr,ncols=1)

for i in range(nr):
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
        try: 
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
        except: 
            print('Bad line in file: ', fn_slug[i])

    print('Finished reading file: ',fn_slug[i])

    tidx = np.argmax(t == tstop) - 1
    #tidx = 1
    ax[i].plot(x, dens[tidx,:], 'r-',\
               x,  vel[tidx,:], 'g-',\
               x, pres[tidx,:], 'b--')
    ax[i].legend(labels=[r'$\rho$',r'$u$',r'$p$'])
    ax[i].set_xlim([0,1])
    ax[i].set_ylim([-3,3])

    ax[i].set_title(sbplt_title[i])

plt.savefig('part2_sol.png',dpi=800)

    
