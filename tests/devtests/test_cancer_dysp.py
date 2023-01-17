"""
Script to explore linking cancer probability to dysplasia

"""

# Imports
import numpy as np
import sciris as sc
import hpvsim as hpv
import hpvsim.utils as hpu
import matplotlib.pyplot as plt
from scipy.stats import lognorm


# Create sim to get baseline prognoses parameters
sim = hpv.Sim(genotypes=[16,18,'hrhpv'])
# sim = hpv.Sim(genotypes=[16,18,31,33,35,51,52,56,58])
sim.initialize()

# Get parameters
ng = sim['n_genotypes']
genotype_pars = sim['genotype_pars']
genotype_map = sim['genotype_map']
cancer_thresh = 0.99


# Shorten duration names
dur_precin = [genotype_pars[genotype_map[g]]['dur_precin'] for g in range(ng)]
dur_dysp = [genotype_pars[genotype_map[g]]['dur_dysp'] for g in range(ng)]
dysp_rate = [genotype_pars[genotype_map[g]]['dysp_rate'] for g in range(ng)]
prog_rate = [genotype_pars[genotype_map[g]]['prog_rate'] for g in range(ng)]
prog_rate_sd = [genotype_pars[genotype_map[g]]['prog_rate_sd'] for g in range(ng)]


def lognorm_params(par1, par2):
    """
    Given the mean and std. dev. of the log-normal distribution, this function
    returns the shape and scale parameters for scipy's parameterization of the
    distribution.
    """
    mean = np.log(par1 ** 2 / np.sqrt(par2 ** 2 + par1 ** 2))  # Computes the mean of the underlying normal distribution
    sigma = np.sqrt(np.log(par2 ** 2 / par1 ** 2 + 1))  # Computes sigma for the underlying normal distribution

    scale = np.exp(mean)
    shape = sigma
    return shape, scale


def logf1(x, k):
    '''
    The concave part of a logistic function, with point of inflexion at 0,0
    and upper asymptote at 1. Accepts 1 parameter which determines the growth rate.
    '''
    return (2 / (1 + np.exp(-k * x))) - 1

def set_font(size=None, font='Libertinus Sans'):
    ''' Set a custom font '''
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return


#%% Functions
def run_calcs():

    # Group genotypes
    genotypes = ['hpv16', 'hpv18', 'hrhpv']
    sim = hpv.Sim(genotypes=genotypes)
    sim.initialize()
    sim['genotype_pars']['hpv16']['prog_rate'] = 0.3
    sim['genotype_pars']['hpv18']['prog_rate'] = 0.4

    # Get parameters
    ng = sim['n_genotypes']
    genotype_map = sim['genotype_map']

    # Get parameters
    genotype_pars = sim['genotype_pars']

    # Shorten duration names
    dur_precin = [genotype_pars[genotype_map[g]]['dur_precin'] for g in range(ng)]
    dur_dysp = [genotype_pars[genotype_map[g]]['dur_dysp'] for g in range(ng)]
    dysp_rate = [genotype_pars[genotype_map[g]]['dysp_rate'] for g in range(ng)]
    prog_rate = [genotype_pars[genotype_map[g]]['prog_rate'] for g in range(ng)]
    cancer_probs = [0.001, 0.002, 0.0005] # Placeholders


    set_font(size=32)
    colors = sc.gridcolors(ng)
    cmap = plt.cm.Oranges([0.25, 0.5, 0.75, 1])
    fig, ax = plt.subplot_mosaic('A;B;C', figsize=(16, 20))


    ####################
    # Make plots
    ####################

    thisx = np.linspace(0.01, 25, 100)
    n_samples = 10

    def cum_cancer_prob(cp,x,dysp): return 1 - np.power(1-(1-np.power(1-cp,dysp*100)),x)

    # Durations and severity of dysplasia
    for gi, gtype in enumerate(genotypes):
        sigma, scale = lognorm_params(dur_dysp[gi]['par1'], dur_dysp[gi]['par2'])
        rv = lognorm(sigma, 0, scale)
        ax['A'].plot(thisx, rv.pdf(thisx), color=colors[gi], lw=2, label=gtype.upper())
        ax['B'].plot(thisx, logf1(thisx, prog_rate[gi]), color=colors[gi], lw=3, label=gtype.upper())
        for smpl in range(n_samples):
            pr = hpu.sample(dist='normal', par1=prog_rate[gi], par2=prog_rate_sd[gi])
            ax['B'].plot(thisx, logf1(thisx, pr), color=colors[gi], lw=1, alpha=0.5, label=gtype.upper())

        cp = cum_cancer_prob(cancer_probs[gi],thisx,logf1(thisx, prog_rate[gi]))
        ax['C'].plot(thisx, cp, color=colors[gi], lw=2, label=gtype.upper())


    ax['A'].set_ylabel("")
    ax['A'].grid()
    ax['A'].set_ylabel("Density")

    ax['B'].set_ylabel("Degree of dysplasia")
    ax['B'].set_ylim([0, 1])
    ax['B'].axhline(y=0.33, ls=':', c='k')
    ax['B'].axhline(y=0.67, ls=':', c='k')
    ax['B'].axhspan(0, 0.33, color=cmap[0], alpha=.4)
    ax['B'].axhspan(0.33, 0.67, color=cmap[1], alpha=.4)
    ax['B'].axhspan(0.67, 1, color=cmap[2], alpha=.4)
    ax['B'].text(-0.3, 0.08, 'CIN1', rotation=90)
    ax['B'].text(-0.3, 0.4, 'CIN2', rotation=90)
    ax['B'].text(-0.3, 0.73, 'CIN3', rotation=90)

    ax['C'].set_xlabel("Duration of dysplasia prior to\nregression/cancer (years)")
    ax['C'].set_ylabel("Probability of cervical\ncancer invasion")
    ax['C'].legend(fontsize=20, frameon=True, loc='best')


    fig.tight_layout()
    plt.savefig(f"AA_cells.png", dpi=100)


#%% Run as a script
if __name__ == '__main__':

    T = sc.tic()

    run_calcs()

    sc.toc(T)
    print('Done.')
