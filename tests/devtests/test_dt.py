import sciris as sc
import pylab as pl
import hpvsim as hpv

# Define the parameters
pars = dict(
    n_agents      = 5e3,       # Population size
    start         = 1980,       # Starting year
    n_years       = 50,         # Number of years to simulate
    rand_seed     = 2,          # Set a non-default seed
    genotypes     = [16, 18],   # Include the two genotypes of greatest general interest
    rel_init_prev = 4.0,
)

# Create the sim
dts = [1/12, 3/12, 4/12, 6/12, 12/12] # use timesteps that produce integral multiples of 12
n_sims = 10
msims = []
for dt in dts:
    sims = []
    for s in range(n_sims):
        pars['rand_seed'] = s
        sim = hpv.Sim(pars, dt=dt, label=f'dt={dt}')
        sims.append(sim)
    msim = hpv.MultiSim(sims)
    msim.run()
    msim.mean()
    msims.append(msim)

merged = hpv.MultiSim.merge(msims, base=True)
merged.plot(['n_alive', 'n_infected', 'hpv_prevalence'], color_by_sim=True)