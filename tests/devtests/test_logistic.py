'''
Tests for single simulations
'''

#%% Imports and settings
import sciris as sc
import numpy as np
import hpvsim as hpv
import pylab as pl
import hpvsim.parameters as hppar
import hpvsim.utils as hpu
import hpvsim as hpv



#%% Run as a script
if __name__ == '__main__':

    # Start timing
    T = sc.tic()


    def set_font(size=None, font='Libertinus Sans'):
        ''' Set a custom font '''
        sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
        sc.options(font=font, fontsize=size)
        return

    set_font(16)

    fig, axes = pl.subplots(3, 1, figsize=(8, 12))
    colors = sc.gridcolors(10)
    dt = 0.01
    t = np.arange(0,30,dt) # Array of years
    x_infl = 7 # Fix point of inflection
    ttc = 15
    tp = 5 / 1e5
    s = 1
    rel_sev=1
    # Try different growth rates:
    # karr = np.linspace(0.001, 1, 6)
    k = 0.4
    for ittc, ttc in enumerate([5, 10, 25]):
        ax = axes[0]
        dysp = hpu.logf3(t * rel_sev, k, x_infl, ttc, s)
        ax.plot(t, dysp, color=colors[ittc], label=f'ttc={ttc}')
        ax.set_title(f'Severity')
        ax.legend()

        ax = axes[1]
        cum_dysp = hppar.compute_severity_integral(t, rel_sev=rel_sev, pars=dict(form='logf3', k=k, x_infl=x_infl, s=s, ttc=ttc))
        cumsum_dysp = np.cumsum(dysp*dt)
        ax.plot(t, cum_dysp, color=colors[ittc], label='integral')
        ax.plot(t, cumsum_dysp, color=colors[ittc], linestyle='--', label='cumsum')
        ax.legend(frameon=False)

        ax = axes[2]
        y1 = hpu.transform_prob(tp, cum_dysp)
        ax.plot(t, y1, color=colors[ittc])
        y2 = hpu.transform_prob(tp, cumsum_dysp)
        ax.plot(t, y2, color=colors[ittc], linestyle='--')

        ax.set_title('Probability of transformation')

    pl.show()

    fig.savefig(f'sevs.png')

    sc.toc(T)
    print('Done.')
