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


# # Create sim to get baseline prognoses parameters
# sim = hpv.Sim(genotypes=[16,18,'hrhpv'])
# # sim = hpv.Sim(genotypes=[16,18,31,33,35,51,52,56,58])
# sim.initialize()
#
# # Get parameters
# ng = sim['n_genotypes']
# genotype_pars = sim['genotype_pars']
# genotype_map = sim['genotype_map']
# cancer_thresh = 0.99
#
#
# # Shorten duration names
# dur_precin = [genotype_pars[genotype_map[g]]['dur_precin'] for g in range(ng)]
# dur_dysp = [genotype_pars[genotype_map[g]]['dur_dysp'] for g in range(ng)]
# dysp_rate = [genotype_pars[genotype_map[g]]['dysp_rate'] for g in range(ng)]
# prog_rate = [genotype_pars[genotype_map[g]]['prog_rate'] for g in range(ng)]
# prog_rate_sd = [genotype_pars[genotype_map[g]]['prog_rate_sd'] for g in range(ng)]


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


def logf2(x, x_infl, k):
    '''
    Logistic function, constrained to pass through 0,0 and with upper asymptote
    at 1. Accepts 2 parameters: growth rate and point of inflexion.
    '''
    l_asymp = -1/(1+np.exp(k*x_infl))
    return l_asymp + 1/( 1 + np.exp(-k*(x-x_infl)))

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

    # Get parameters
    ng = sim['n_genotypes']
    genotype_map = sim['genotype_map']

    # Get parameters
    genotype_pars = sim['genotype_pars']
    genotype_pars['hpv16']['prog_rate'] = 0.2
    genotype_pars['hpv18']['prog_rate'] = 0.3

    genotype_pars['hpv16']['dysp_rate'] = 0.5
    genotype_pars['hpv18']['dysp_rate'] = 0.4
    genotype_pars['hrhpv']['dysp_rate'] = 0.1

    genotype_pars['hpv16']['dysp_infl'] = 6
    genotype_pars['hpv18']['dysp_infl'] = 6
    genotype_pars['hrhpv']['dysp_infl'] = 6

    genotype_pars['hpv16']['dur_precin']['par1'] = 2
    genotype_pars['hpv18']['dur_precin']['par1'] = 1.5
    genotype_pars['hrhpv']['dur_precin']['par1'] = 2.5
    genotype_pars['hpv16']['dur_precin']['par2'] = 2
    genotype_pars['hpv18']['dur_precin']['par2'] = 1.5
    genotype_pars['hrhpv']['dur_precin']['par2'] = 3


    # Shorten duration names
    dur_prod = [genotype_pars[genotype_map[g]]['dur_precin'] for g in range(ng)]
    dur_trans = [genotype_pars[genotype_map[g]]['dur_dysp'] for g in range(ng)]
    trans_rate = [genotype_pars[genotype_map[g]]['dysp_rate'] for g in range(ng)]
    trans_infl = [genotype_pars[genotype_map[g]]['dysp_infl'] for g in range(ng)]
    prog_rate = [genotype_pars[genotype_map[g]]['prog_rate'] for g in range(ng)]
    prog_rate_sd = [genotype_pars[genotype_map[g]]['prog_rate_sd'] for g in range(ng)]
    cancer_probs = [0.0005, 0.0005, 0.0005] # Placeholders


    set_font(size=20)
    colors = sc.gridcolors(ng)
    cmap = plt.cm.Oranges([0.25, 0.5, 0.75, 1])
    fig, ax = plt.subplot_mosaic('AB;CD;EF', figsize=(16, 20))

    ####################
    # Panel A and C
    ####################

    x = np.linspace(0.01, 10, 200)  # Make an array of durations 0-3 years
    glabels = ['HPV16', 'HPV18', 'HRHPV']
    dysp_shares = []
    gtypes = []
    igi = 0.01  # Define the integration interval
    longx = sc.inclusiverange(0.01, 20, igi)  # Initialize a LONG array of years

    # Loop over genotypes, plot each one
    for gi, gtype in enumerate(genotypes):
        sigma, scale = lognorm_params(dur_prod[gi]['par1'], dur_prod[gi]['par2'])
        rv = lognorm(sigma, 0, scale)
        aa = np.diff(rv.cdf(longx))  # Calculate the probability that a woman will have a pre-dysplasia duration in any of the subintervals of time spanning 0-25 years
        bb = logf2(longx, trans_infl[gi], trans_rate[gi])[1:]  # Calculate the probablity of her developing dysplasia for a given duration
        dysp_shares.append(np.dot(aa,bb))  # Convolve the two above calculations to determine the probability of her developing dysplasia overall
        gtypes.append(gtype)  # Store genotype names for labeling
        ax['A'].plot(x, rv.pdf(x), color=colors[gi], lw=2, label=glabels[gi])
        ax['C'].plot(x, logf2(x, trans_infl[gi], trans_rate[gi]), color=colors[gi], lw=2, label=gtype.upper())

    bottom = np.zeros(ng)
    ax['E'].bar(np.arange(1, ng + 1), 1-np.array(dysp_shares), color='grey', bottom=bottom, label='Productive')
    ax['E'].bar(np.arange(1, ng + 1), np.array(dysp_shares), color=cmap[0], bottom=1-np.array(dysp_shares), label='Transforming')
    ax['E'].set_xticks(np.arange(1, ng + 1))
    ax['E'].set_xticklabels(glabels)
    ax['E'].set_ylabel("")
    ax['E'].set_ylabel("Distribution of infection outcomes")

    # Axis labeling and other settings
    ax['C'].set_xlabel("Duration of productive infection (years)")
    for axn in ['A', 'C']:
        ax[axn].set_ylabel("")
        ax[axn].grid()

    ax['A'].set_ylabel("Density")
    ax['C'].set_ylabel("Probability of transformation")
    ax['A'].set_xlabel("Duration of productive infection prior to\nclearance or transformation (years)")

    ax['A'].legend(fontsize=20, frameon=False)
    ax['E'].legend(fontsize=20, frameon=True)


    ####################
    # Make plots
    ####################

    thisx = np.linspace(0.01, 25, 100)
    n_samples = 10

    def cum_cancer_prob(cp,x,dysp): return 1 - np.power(1-(1-np.power(1-cp,dysp*100)),x)

    # Durations and severity of dysplasia
    for gi, gtype in enumerate(genotypes):
        sigma, scale = lognorm_params(dur_trans[gi]['par1'], dur_trans[gi]['par2'])
        rv = lognorm(sigma, 0, scale)
        ax['B'].plot(thisx, rv.pdf(thisx), color=colors[gi], lw=2, label=gtype.upper())
        ax['D'].plot(thisx, logf1(thisx, prog_rate[gi]), color=colors[gi], lw=3, label=gtype.upper())
        for smpl in range(n_samples):
            pr = hpu.sample(dist='normal', par1=prog_rate[gi], par2=prog_rate_sd[gi])
            ax['D'].plot(thisx, logf1(thisx, pr), color=colors[gi], lw=1, alpha=0.5, label=gtype.upper())

        cp = cum_cancer_prob(cancer_probs[gi],thisx,logf1(thisx, prog_rate[gi]))
        ax['F'].plot(thisx, cp, color=colors[gi], lw=2, label=gtype.upper())


    ax['B'].set_ylabel("")
    ax['B'].grid()
    ax['B'].set_ylabel("Density")

    ax['B'].set_xlabel("Duration of transforming infection prior to\nregression/cancer (years)")
    ax['D'].set_xlabel("Duration of transforming infection (years)")

    ax['D'].set_ylabel("Degree of transformation")
    ax['D'].set_ylim([0, 1])
    ax['D'].axhline(y=0.5, ls=':', c='k')
    ax['D'].axhspan(0, 0.5, color=cmap[1], alpha=.4)
    ax['D'].axhspan(0.5, 1, color=cmap[2], alpha=.4)
    ax['D'].text(-0.3, 0.22, 'CIN2', rotation=90)
    ax['D'].text(-0.3, 0.72, 'CIN3', rotation=90)

    ax['F'].set_xlabel("Duration of transforming infection (years)")
    ax['F'].set_ylabel("Probability of cervical\ncancer invasion")
    ax['F'].legend(fontsize=20, frameon=True, loc='best')
    # ax['E'].set_axis_off()


    fig.tight_layout()
    plt.savefig(f"AA_cells.png", dpi=100)


#%% Run as a script
if __name__ == '__main__':

    T = sc.tic()

    run_calcs()

    sc.toc(T)
    print('Done.')
