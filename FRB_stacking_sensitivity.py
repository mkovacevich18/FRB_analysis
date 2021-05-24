### This script serves to mainly run the background and signal trials and save them as a nested dictionary, it is optimized to run on the npx cluster ###

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import argparse
import histlite as hl
import csky as cy

ana_dir = cy.utils.ensure_dir('/data/user/mkovacevich/FRB_analysis/cascades_ana')
repo = cy.selections.Repository()
ana = cy.analysis.Analysis(repo, cy.selections.MESEDataSpecs.mesc_7yr, dir=ana_dir)

trials_dir = cy.utils.ensure_dir('/data/user/mkovacevich/FRB_analysis/trials')
sig_dir = cy.utils.ensure_dir('{}/sig'.format(trials_dir))
bg_dir = cy.utils.ensure_dir('{}/bg'.format(trials_dir))

parser = argparse.ArgumentParser(description='Process sensitivities for FRB catalog over livetime of MESC 7yr dataset')
parser.add_argument('--gamma',type=float,help='Spectral indice for E')
parser.add_argument('--dt', type=float, help='time-window in units of seconds')
parser.add_argument('--seed',type=float,help ='keep track of the total number of jobs submitted')
parser.add_argument('--choose',type=float,help='decide between bg (1), sig (2) and flux calculations (3)')
parser.add_argument('--bg_trials',type=float,help='the number of bg trials to run')
parser.add_argument('--sig_trials',type=float,help='the number of signal trials to run')
parser.add_argument('--cpus',type=float,help='choose the number of cpus to run for a job')
args = parser.parse_args()

cy.CONF['ana'] = ana
cy.CONF['mp_cpus'] = args.cpus

#MJD, RA and Dec of each FRB
FRB_mjd_time = [57760.83732799769, 57547.161818148146, 57488.35670925926, 57464.37542280093, 57389.35323349537, 57386.67762181713, 57362.262416412035, 57241.74578506944, 57183.22707634259, 57130.18688260417, 57068.86228835648, 56791.718183564815, 56600.752907407405, 56502.376286921295, 56471.165279837965, 56469.62221957176, 56202.5481300463, 55953.341223668984, 55745.79144221065, 55738.89811894676, 55704.629394675925, 55612.080417777775]

FRB_ra_rad = [2.979102500229121, 1.992642407001926, 2.2750366799746087, 2.0673424989872835, 5.928883469024737, 2.534392606820966, 5.067563483165536, 5.945115031068284, 2.811899957888065, 1.9050268785518107, 4.788485335771643, 5.9082885838512045, 1.7634806762150708, 3.5838590860451562, 2.369458992507502, 4.306948995146407, 4.776966162708479, 6.087184832180623, 6.1559508047091995, 5.514018705825685, 5.695009349257497, 5.910732044803997]

FRB_dec_rad = [-0.08691739674931762, -0.6841690667817772, 0.10611601852125524, -0.4954989746411902, -0.5204571829447091, -0.0445058959258554, -0.06754424205218056, -0.9585348251952858, -0.6965609044709369, -0.33143802495372315, -0.054105206811824215, -0.2040289895581371, -0.8852309966115239, -0.06981317007977318, 0.060039326268604934, -0.11414453308042916, -1.4800392056911915, -0.3066543495754037, -0.01972222054753592, -0.755029434412747, 0.00027925268031909274, -0.20245819323134223]

#Length of the time window after a burst in which to search for neutrinos
FRB_time_window = np.ones_like(FRB_ra_rad)*args.dt/86400.

#Defining the number of signal and background trials to run (and how many signal events to inject/interpolate over)
if args.dt > 1e3:
    n_bg_trials = args.bg_trials/1e2
    
    if args.dt > 1e5:
        n_sigs = np.r_[2:10:1, 40:100.1:2]
    elif args.dt <= 1e5:
        n_sigs = np.r_[2:10:1, 10:30.1:2]
        
elif args.dt <= 1e3:
    n_bg_trials = args.bg_trials
    n_sigs = np.r_[2:10:1, 10:30.1:2]
        
n_sig_trials = args.sig_trials


##### Defining functions to compute background trials, signal trials and flux ######    

def do_background_trials(N=n_bg_trials):
    src = cy.sources(FRB_ra_rad, FRB_dec_rad, mjd = FRB_mjd_time, sigma_t = np.zeros_like(FRB_ra_rad), t_100 = FRB_time_window)
    conf = {'extended':True, 'space':"ps",'time':"transient",'sig':"transient",'flux': cy.hyp.PowerLawFlux(args.gamma)}
    tr = cy.get_trial_runner(conf, src = src, ana=ana)
    # run trials
    trials = tr.get_many_fits(N,logging=False)
    # save to disk
    dir = cy.utils.ensure_dir('{}/dt/{}'.format(bg_dir, args.dt))
    filename = '{}/bg_trials_seed_{}.npy'.format(dir, args.seed)
    print('->', filename)
    np.save(filename, trials.as_array)
    
#n_sigs = np.r_[2:10:1, 40:100.1:2] extend to 100.1 instead of 30.1 if analyzing time windows greater than 1e5 seconds
def do_signal_trials(n_sig, N=n_sig_trials):
    src = cy.sources(FRB_ra_rad, FRB_dec_rad, mjd = FRB_mjd_time, sigma_t = np.zeros_like(FRB_ra_rad), t_100 = FRB_time_window)
    conf = {'extended':True, 'space':"ps",'time':"transient",'sig':"transient",'flux': cy.hyp.PowerLawFlux(args.gamma)} 
    tr = cy.get_trial_runner(conf, src = src, ana=ana)
    trials = tr.get_many_fits(N, n_sig, poisson = True, logging=False)
    dir = cy.utils.ensure_dir('{}/gamma/{}/dt/{}/n_sig/{}'.format(sig_dir, args.gamma, args.dt, n_sig))
    filename = '{}/sig_trials_{}.npy'.format(dir, args.seed)
    print('->', filename)
    np.save(filename, trials.as_array)
    
def ndarray_to_TSD(trials):
    return cy.dists.TSD(cy.utils.Arrays(trials))

def find_n_sig(beta=0.9, nsigma=None):
    sig_trials = cy.bk.get_best(sig, 'gamma', args.gamma, 'dt', args.dt, 'n_sig')
    b = cy.bk.get_best(bg, 'gamma', args.gamma, 'dt', args.dt)
    src = cy.sources(FRB_ra_rad, FRB_dec_rad, mjd = FRB_mjd_time, sigma_t = np.zeros_like(FRB_ra_rad), t_100 = FRB_time_window)
    conf = {'extended':True, 'space':"ps",'time':"transient",'sig':"transient",'flux': cy.hyp.PowerLawFlux(args.gamma)}
    tr = cy.get_trial_runner(conf, src = src, ana=ana)
    if nsigma is not None:
        ts = b.isf_nsigma(nsigma)
    else:
        ts = b.median()
    trials = {0: b.trials}
    trials.update(sig_trials)
    result = tr.find_n_sig(ts, beta, max_batch_size=0, logging=False, trials=trials, n_bootstrap=1)
    return tr.to_E2dNdE(result, E0=1e5)

#Computing background or signal trials and saving them to /data/user/mkovacevich
if args.choose == 1:
    print('Computing bg trials w/ the time window ' + str(args.dt))
    do_background_trials(N=n_bg_trials)
    
elif args.choose == 2:
    print('Computing signal trials for gamma = ' + str(args.gamma) + 'and w/ the time window ' + str(args.dt))
    for n_sig in n_sigs:
        do_signal_trials(n_sig, N=n_sig_trials)
        
### This option is not used in this script, it is not computationally intensive to calculate n_sig/sensitivities/disc potentials once the trials have been performed and saved. Instead we use the notebook ("FRB_stacking_sensitivity_flux.ipynb") which utilizes this same layout but uses the cobalts rather than the npx cluster
elif args.choose == 3:
    bg = cy.bk.get_all('{}/'.format(bg_dir),'bg*npy', merge=np.concatenate,post_convert=ndarray_to_TSD)
    sig = cy.bk.get_all('{}/'.format(sig_dir),'sig*npy', merge=np.concatenate, post_convert=cy.utils.Arrays)
    fluxs_sens = []
    fluxs_sens = find_n_sig(beta=0.5)
    print('sens = ' + str(fluxs_sens))
    np.save('sens_gamma_'+str(args.gamma)+'dt'+str(args.dt),fluxs_sens,allow_pickle=True)
    
    fluxs_3sig_disc,fluxs_4sig_disc, fluxs_5sig_disc = [],[],[]
    fluxs_3sig_disc = find_n_sig(beta=0.9,nsigma=3)
    print('3sig = ' + str(fluxs_3sig_disc))
    np.save('3sig_gamma_'+str(args.gamma)+'dt'+str(args.dt),fluxs_3sig_disc,allow_pickle=True)
    
    fluxs_4sig_disc = find_n_sig(beta=0.9,nsigma=4)
    np.save('4sig_gamma_'+str(args.gamma)+'dt'+str(args.dt),fluxs_4sig_disc,allow_pickle=True)
    print('4sig = ' + str(fluxs_4sig_disc))
    
    fluxs_5sig_disc = find_n_sig(beta = 0.9, nsigma=5)
    print('5sig = ' + str(fluxs_5sig_disc))
    np.save('5sig_gamma_'+str(args.gamma)+'dt'+str(args.dt),fluxs_5sig_disc,allow_pickle=True)

    
