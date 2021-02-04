'''
Programmer: Mike Kovacevich
Last edited: 1/11/21
E-mail: mgk56@drexel.edu

Differential stacking sensitivity - objective is to compute sensitivity vs. energy for different time windows and gamma 
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import argparse
import histlite as hl
import csky as cy
import random

cy.plotting.mrichman_mpl()

timer = cy.timing.Timer()
time = timer.time

parser = argparse.ArgumentParser(description='Process sensitivities for FRB catalog over livetime of MESC 7yr dataset')
parser.add_argument('--gamma',type=float,help='Spectral indice for E')
parser.add_argument('--dt', type=float, help='time-window in units of seconds')
#parser.add_argument('--sigma',type=float, help = 'which discovery potential (3,4 or 5 sigma @ 90%) to look at')
parser.add_argument('--cpus',type=float,help='choose the number of cpus to run for a job')
args = parser.parse_args()

gamma = args.gamma
dt = args.dt
n_cpus = args.cpus

FRB_mjd_time = [57760.83732799769, 57547.161818148146, 57488.35670925926, 57464.37542280093, 57389.35323349537, 57386.67762181713, 57362.262416412035, 57241.74578506944, 57183.22707634259, 57130.18688260417, 57068.86228835648, 56791.718183564815, 56600.752907407405, 56502.376286921295, 56471.165279837965, 56469.62221957176, 56202.5481300463, 55953.341223668984, 55745.79144221065, 55738.89811894676, 55704.629394675925, 55612.080417777775]

FRB_ra_rad = [2.979102500229121, 1.992642407001926, 2.2750366799746087, 2.0673424989872835, 5.928883469024737, 2.534392606820966, 5.067563483165536, 5.945115031068284, 2.811899957888065, 1.9050268785518107, 4.788485335771643, 5.9082885838512045, 1.7634806762150708, 3.5838590860451562, 2.369458992507502, 4.306948995146407, 4.776966162708479, 6.087184832180623, 6.1559508047091995, 5.514018705825685, 5.695009349257497, 5.910732044803997]

FRB_dec_rad = [-0.08691739674931762, -0.6841690667817772, 0.10611601852125524, -0.4954989746411902, -0.5204571829447091, -0.0445058959258554, -0.06754424205218056, -0.9585348251952858, -0.6965609044709369, -0.33143802495372315, -0.054105206811824215, -0.2040289895581371, -0.8852309966115239, -0.06981317007977318, 0.060039326268604934, -0.11414453308042916, -1.4800392056911915, -0.3066543495754037, -0.01972222054753592, -0.755029434412747, 0.00027925268031909274, -0.20245819323134223]

FRB_time_window = np.ones_like(FRB_ra_rad)*dt/86400.

#Building/loading MESC data from analysis directory
ana_dir = cy.utils.ensure_dir('/data/user/mkovacevich/FRB_analysis/cascades_ana')
repo = cy.selections.Repository()
ana = cy.analysis.Analysis(repo, cy.selections.MESEDataSpecs.mesc_7yr, dir=ana_dir)

cy.CONF['ana'] = ana
cy.CONF['mp_cpus'] = n_cpus

src = cy.sources(FRB_ra_rad, FRB_dec_rad, mjd = FRB_mjd_time, sigma_t = np.zeros_like(FRB_ra_rad), t_100 = FRB_time_window)
conf = {'extended':True, 'space':"ps",'time':"transient",'sig':"transient"} #,'flux': cy.hyp.PowerLawFlux(gamma)}

with time('background estimation'):
    allE_tr = cy.get_trial_runner(conf, ana=ana, src=src)
    bg = cy.dists.TSD(allE_tr.get_many_fits(10000, mp_cpus=n_cpus))
    
# you could also use np.logspace()
Ebins = 10**np.r_[3:7.1:.25]
trs = [
    cy.get_trial_runner(conf, 
        ana=ana, src=src,
        flux=cy.hyp.PowerLawFlux(gamma, energy_range=(Emin, Emax)))
    for (Emin, Emax) in zip(Ebins[:-1], Ebins[1:])
]

senss = []

with time('differential sensitivity'):
    for (Emin, Emax, tr) in zip(Ebins[:-1], Ebins[1:], trs):
        print(f'log10(E) bin: {np.log10(Emin), np.log10(Emax)}')
        sens = tr.find_n_sig(
            0, .9, batch_size=500, max_batch_size=500,
            seed=random.randint(0, 1e3), mp_cpus=n_cpus)
        # N.B. must set E0 to a value within the energy_range
        # of the flux for this trial runner
        sens['E2dNdE'] = tr.to_E2dNdE(sens, E0=Emin)
        senss.append(sens)
        
np.save('/data/user/mkovacevich/FRB_analysis/DiffSens_gamma_'+str(gamma)+'_dt_'+str(dt)+'.npy', senss, allow_pickle = True)
'''
fluxs = [s['E2dNdE']/1e3 for s in senss] # convert GeV->TeV
diffsens = hl.Hist(Ebins, fluxs)
#diffsens

fig, ax = plt.subplots()
hl.plot1d(ax, diffsens)
ax.loglog()
ax.set(
    ylim=1e-3,
    xlabel=r'$E$ [GeV]',
    ylabel=r'$E^2\cdot dN/dE~~[\text{TeV}\,\text{cm}^2\,\text{}]$',
    title=r'Differential Sensitivity (dt = ' + str(dt) +' s and gamma = ' + str(gamma) +')'
)
ax.grid()
fig.savefig('/data/user/mkovacevich/FRB_analysis/DifferentialSensitivity_gamma_'+str(gamma)+'_dt_'+str(dt)+'.png')
'''