import hpvsim as hpv
import matplotlib.pyplot as plt

# Define the parameters
pars = dict(
    n_agents      = 5e3,       # Population size
    start         = 1980,       # Starting year
    n_years       = 10,         # Number of years to simulate
    rand_seed     = 2,          # Set a non-default seed
    genotypes     = [16, 18],   # Include the two genotypes of greatest general interest
    rel_init_prev = 4.0,
)

# Create the sim
dts = [0.625, 0.125, 0.25, 0.5, 1]

n_sims = 10
msims = []
for dt in dts:
    sims = []
    for s in range(n_sims):
        pars['rand_seed'] = s
        sim = hpv.Sim(pars, dt=dt, label=f'dt={dt}')
        #sim.run()
        sims.append(sim)
    msim = hpv.MultiSim(sims)
    msim.run()
    msim.median()
    msims.append(msim)

merged = hpv.MultiSim.merge(msims, base=True)

sim_base = merged.sims[-1]

for sim in merged.sims:
    plt.plot(sim_base.results['year'], sim.results['infections'].values/sim_base.results['infections'].values)

plt.xlabel('year')
plt.ylabel('ratio dt/dt=1 infections')
plt.show()

