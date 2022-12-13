'''
Defines classes and methods for calculating immunity
'''

import numpy as np
import sciris as sc
from collections.abc import Iterable
from . import utils as hpu
from . import defaults as hpd
from . import parameters as hppar
from . import interventions as hpi


# %% Immunity methods

def init_immunity(sim, create=True):
    ''' Initialize immunity matrices with all genotypes and vaccines in the sim'''

    # Pull out all of the circulating genotypes for cross-immunity
    ng = sim['n_genotypes']

    # Pull out all the vaccination interventions
    vx_intvs = [x for x in sim['interventions'] if isinstance(x, hpi.BaseVaccination)]
    nv = len(vx_intvs)

    # Dimension for immunity matrix
    ndim = ng + nv

    # If immunity values have been provided, process them
    if sim['immunity'] is None or create:

        # Precompute waning - same for all genotypes
        if sim['use_waning']:
            imm_decay = sc.dcp(sim['imm_decay'])
            if 'half_life' in imm_decay.keys():
                imm_decay['half_life'] /= sim['dt']
            sim['imm_kin'] = precompute_waning(t=sim.tvec, pars=imm_decay)

        sim['immunity_map'] = dict()
        # Firstly, initialize immunity matrix with defaults. These are then overwitten with specific values below
        immunity = np.ones((ng, ng), dtype=hpd.default_float)  # Fill with defaults

        # Next, overwrite these defaults with any known immunity values about specific genotypes
        default_cross_immunity = hppar.get_cross_immunity(cross_imm_med=sim['cross_imm_med'], cross_imm_high=sim['cross_imm_high'])
        for i in range(ng):
            sim['immunity_map'][i] = 'infection'
            label_i = sim['genotype_map'][i]
            for j in range(ng):
                label_j = sim['genotype_map'][j]
                if label_i in default_cross_immunity and label_j in default_cross_immunity:
                    immunity[j][i] = default_cross_immunity[label_j][label_i]

        imm_source = ng
        for vi,vx_intv in enumerate(vx_intvs):
            genotype_pars_df = vx_intv.product.genotype_pars[vx_intv.product.genotype_pars.genotype.isin(sim['genotype_map'].values())] # TODO fix this
            vacc_mapping = [genotype_pars_df[genotype_pars_df.genotype==gtype].rel_imm.values[0] for gtype in sim['genotype_map'].values()]
            vacc_mapping += [1]*(vi+1) # Add on some ones to pad out the matrix
            vacc_mapping = np.reshape(vacc_mapping, (len(immunity)+1, 1)).astype(hpd.default_float) # Reshape
            immunity = np.hstack((immunity, vacc_mapping[0:len(immunity),]))
            immunity = np.vstack((immunity, np.transpose(vacc_mapping)))
            vx_intv.product.imm_source = imm_source
            imm_source += 1

        sim['immunity'] = immunity

    sim['immunity'] = sim['immunity'].astype('float32')
    sim['n_imm_sources'] = ndim

    return


def update_peak_immunity(people, inds, imm_pars, imm_source, offset=None, infection=True):
    '''
        Update immunity level

        This function updates the immunity for individuals when an infection or vaccination occurs.
            - individuals that are infected and already have immunity from a previous vaccination/infection have their immunity level;
            - individuals without prior immunity are assigned an initial level drawn from a distribution. This level
                depends on whether the immunity is from a natural infection or from a vaccination (and if so, on the type of vaccine).

        Args:
            people: A people object
            inds: Array of people indices
            imm_pars: Parameters from which to draw values for quantities like ['imm_init'] - either
                        sim pars (for natural immunity) or vaccine pars
            imm_source: index of either genotype or vaccine where immunity is coming from

        Returns: None
    '''

    if infection:
        # Determine whether individual seroconverts based upon genotype
        genotype_label = imm_pars['genotype_map'][imm_source]
        genotype_pars = imm_pars['genotype_pars'][genotype_label]
        seroconvert_probs = np.full(len(inds), fill_value=genotype_pars.sero_prob)
        is_seroconvert = hpu.binomial_arr(seroconvert_probs)

        # Extract parameters and indices
        has_imm = people.imm[imm_source, inds] > 0
        no_prior_imm_inds = inds[~has_imm]
        prior_imm_inds = inds[has_imm]

        if len(prior_imm_inds):
            boost = genotype_pars['imm_boost']
            people.peak_imm[imm_source, prior_imm_inds] *= is_seroconvert[has_imm] * boost

        if len(no_prior_imm_inds):
            people.peak_imm[imm_source, no_prior_imm_inds] = is_seroconvert[~has_imm] * hpu.sample(**imm_pars['imm_init'], size=len(no_prior_imm_inds))

    else:
        # Vaccination by dose
        dose1_inds = inds[people.doses[inds]==1] # First doses
        dose2_inds = inds[people.doses[inds]==2] # Second doses
        dose3_inds = inds[people.doses[inds]==3] # Third doses
        if imm_pars['doses']>1:
            imm_pars['imm_boost'] = sc.promotetolist(imm_pars['imm_boost'])
        if len(dose1_inds)>0: # Initialize immunity for newly vaccinated people
            people.peak_imm[imm_source, dose1_inds] = hpu.sample(**imm_pars['imm_init'], size=len(dose1_inds))
        if len(dose2_inds) > 0: # Boost immunity for people receiving 2nd dose...
            people.peak_imm[imm_source, dose2_inds] *= imm_pars['imm_boost'][0]
        if len(dose3_inds) > 0:
            people.peak_imm[imm_source, dose3_inds] *= imm_pars['imm_boost'][1]

    base_t = people.t + offset if offset is not None else people.t
    people.t_imm_event[imm_source, inds] = base_t

    return


def check_immunity(people):
    '''
    Calculate people's immunity on this timestep from prior infections.
    As an example, suppose HPV16 and 18 are in the sim, and the cross-immunity matrix is:
        pars['immunity'] = np.array([[1., 0.5],
                                     [0.3, 1.]])
    This indicates that people who've had HPV18 have 50% protection against getting 16, and
    people who've had 16 have 30% protection against getting 18.
    Now suppose we have 3 people, whose immunity levels are
        people.imm = np.array([[0.9, 0.0, 0.0],
                               [0.0, 0.7, 0.0]])
    This indicates that person 1 has a prior HPV16 infection, person 2 has a prior HPV18
    infection, and person 3 has no history of infection.

    In this function, we take the dot product of pars['immunity'] and people.imm to get:
        people.sus_imm = np.array([[0.9 , 0.35, 0.  ],
                                   [0.27, 0.7 , 0.  ]])
    This indicates that the person 1 has high protection against reinfection with HPV16, and
    some (30% of 90%) protection against infection with HPV18, and so on.

    '''
    immunity = people.pars['immunity'] # cross-immunity/own-immunity scalars to be applied to immunity level
    sus_imm = np.dot(immunity,people.imm) # Dot product gives immunity to all genotypes
    people.sus_imm[:] = np.minimum(sus_imm, np.ones_like(sus_imm)) # Don't let this be above 1
    return




#%% Methods for computing waning

def precompute_waning(t, pars=None):
    '''
    Process functional form and parameters into values:

        - 'exp_decay'   : exponential decay. Parameters should be init_val and half_life (half_life can be None/nan)
        - 'linear_decay': linear decay

    A custom function can also be supplied.

    Args:
        length (float): length of array to return, i.e., for how long waning is calculated
        pars (dict): passed to individual immunity functions

    Returns:
        array of length 'length' of values
    '''

    pars = sc.dcp(pars)
    form = pars.pop('form')
    choices = [
        'exp_decay',
    ]

    # Process inputs
    if form is None:
        output = np.ones(len(t), dtype=hpd.default_float)

    elif form == 'exp_decay':
        if pars['half_life'] is None: pars['half_life'] = np.nan
        output = exp_decay(t, **pars)

    elif callable(form):
        output = form(t, **pars)

    else:
        errormsg = f'The selected functional form "{form}" is not implemented; choices are: {sc.strjoin(choices)}'
        raise NotImplementedError(errormsg)

    return output


def exp_decay(t, init_val, half_life):
    '''
    Returns an array of length t with values for the immunity at each time step after recovery
    '''
    decay_rate = np.log(2) / half_life if ~np.isnan(half_life) else 0.
    result = init_val * np.exp(-decay_rate * t, dtype=hpd.default_float)
    return result


def linear_decay(length, init_val, slope):
    ''' Calculate linear decay '''
    result = -slope*np.ones(length)
    result[0] = init_val
    return result


def linear_growth(length, slope):
    ''' Calculate linear growth '''
    return slope*np.ones(length)
