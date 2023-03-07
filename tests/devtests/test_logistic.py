'''
Tests for single simulations
'''

#%% Imports and settings
import sciris as sc
import numpy as np
import hpvsim as hpv
import pylab as pl
import hpvsim.parameters as hppar

def logf2(x, x_infl, k):
    l_asymp = -1/(1+np.exp(k*x_infl))
    return l_asymp + 1/( 1 + np.exp(-k*(x-x_infl)))

def get_asymptotes(ed50, k, ttc, s):
    term1 = (1 + np.exp(k*(x_infl-ttc)))**s
    term2 = (1 + np.exp(k*x_infl))**s
    u_asymp_num = term1*(1-term2)
    u_asymp_denom = term1 - term2
    u_asymp = u_asymp_num / u_asymp_denom
    l_asymp = term1 / (term1 - term2)
    return l_asymp, u_asymp

def transform_prob(tp,dysp):
    return 1-np.power(1-tp, dysp*100)

def cum_transform_prob(tp, t, dysp):
    dd = np.diff(dysp)
    n = len(t)
    result = [1 - np.product([((1 - tp) ** (100 * dd[i])) ** (j - i) for i in range(j)]) for j in range(n)]
    return result

def logf3(x, k, x_infl, ttc=25, s=1):
    l_asymp, u_asymp = get_asymptotes(k, x_infl, ttc, s)
    return np.minimum(1, l_asymp + (u_asymp-l_asymp)/(1+np.exp(k*(x_infl-x)))**s)


def invlogf3(y, k, x_infl, ttc=25, s=1):
    l_asymp, u_asymp = get_asymptotes(k, x_infl, ttc, s)
    part1 = np.log((u_asymp-l_asymp)/(y-l_asymp))/s
    part2 = np.log(np.exp(part1)-1)
    final = 1/k * (k*x_infl - part2)
    return final

def intlogf3(upper, k, x_infl, ttc=25, s=1):
    l_asymp, u_asymp = get_asymptotes(k, x_infl, ttc, s)
    if rel_sev is not None: upper = rel_sev * upper
    val_at_0    = 1/k* ((u_asymp-l_asymp)*np.log(np.exp(k*x_infl)+1))
    val_at_lim  = 1/k* ((u_asymp-l_asymp)*np.log(np.exp(k*(x_infl-upper))+1)) + u_asymp*upper
    return val_at_lim-val_at_0


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

    fig, ax = pl.subplots(1, 1, figsize=(8, 8))
    colors = sc.gridcolors(10)
    celltp = 0.00025
    t = np.arange(0,30,1) # Array of years
    x_infl = 10 # Fix point of inflection
    k = 0.3
    pars = {'form':'logf3', 'k':k, 'x_infl':x_infl}

    # Try different growth rates:
    rel_sev = 1
    y = hppar.compute_severity(t, rel_sev=rel_sev, pars=pars)
    yint = intlogf3(t*rel_sev, k, x_infl, 30, 1)
    tp1 = transform_prob(celltp, y)
    tp2 = transform_prob(celltp, yint)
    tp3 = cum_transform_prob(celltp, t, y)
    tp4 = transform_prob(celltp, t/100)
    ax.plot(t, tp1, color=colors[0], label=f'1-(1-celltp)**(100*dysp)\n with dysp=logf3(dur, {k=},{x_infl=})')
    ax.plot(t, tp2, color=colors[1], label=f'1-(1-celltp)**(100*dysp)\n with dysp=\int_0^dur logf3(u, {k=},{x_infl=})du')
    ax.plot(t, tp3, color=colors[3], label=f'From paper')
    ax.plot(t, tp4, color=colors[4], label=f'1-(1-celltp)**t')
    ax.set_xlabel('Duration of infection')
    ax.set_title(f'Assigned transformation probability with {celltp=}')

    ax.legend(frameon=False)
    pl.show()
    fig.savefig(f'tps.png')

    fig, ax = pl.subplots(1, 1, figsize=(8, 8))
    ax.plot(t, y, color=colors[0], label=f'dysp=logf3(dur, {k=},{x_infl=})')
    ax.set_xlabel('Duration of infection')
    ax.set_title(f'Severity')
    fig.savefig(f'sevs.png')

    sc.toc(T)
    print('Done.')
