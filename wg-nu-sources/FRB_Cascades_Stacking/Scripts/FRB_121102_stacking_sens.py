import numpy as np
import argparse
import histlite as hl
import csky as cy

#Building/loading MESC data from analysis directory
ana_dir = cy.utils.ensure_dir('/data/user/mkovacevich/FRB_analysis/cascades_ana')
repo = cy.selections.Repository()
ana = cy.analysis.Analysis(repo, cy.selections.MESEDataSpecs.mesc_7yr, dir=ana_dir)

#Directories to store the trials
trials_dir = cy.utils.ensure_dir('/data/user/mkovacevich/FRB_analysis/FRB_121102_trials')
sig_dir = cy.utils.ensure_dir('{}/sig'.format(trials_dir))
bg_dir = cy.utils.ensure_dir('{}/bg'.format(trials_dir))

parser = argparse.ArgumentParser(description='Process sensitivities for FRB catalog over livetime of MESC 7yr dataset')
parser.add_argument('--gamma',type=float,help='Spectral indice for E')
parser.add_argument('--dt', type=float, help='time-window in units of seconds')
#parser.add_argument('--sigma',type=float, help = 'which discovery potential (3,4 or 5 sigma @ 90%) to look at')
parser.add_argument('--seed',type=float,help ='keep track of the total number of jobs submitted')
parser.add_argument('--choose',type=float,help='decide between bg (1), sig (2) and flux calculations (3)')
parser.add_argument('--bg_trials',type=float,help='the number of bg trials to run')
parser.add_argument('--sig_trials',type=float,help='the number of signal trials to run')
parser.add_argument('--cpus',type=float,help='choose the number of cpus to run for a job')
args = parser.parse_args()

cy.CONF['ana'] = ana
cy.CONF['mp_cpus'] = args.cpus

FRB_mjd_time = [57623.74402686, 57364.20463266, 57345.4624131, 57345.4575953, 57345.45248793, 57345.44769125, 57339.35604601, 57175.74828727, 57175.74762485, 57175.74566583, 57175.74351039, 57175.74283934, 57175.7425767, 57175.69972782, 57175.69314323, 57159.74422362, 57159.73760083, 56233.27492181, 57633.67986367, 57633.69515937, 57638.49937435, 57643.45730263, 57645.42958602, 57646.4660065, 57648.4369149 , 57649.45175697]

FRB_ra_rad = [1.44852874, 1.44847784, 1.44847784, 1.44847784, 1.44847784, 1.44847784, 1.44847784, 1.44825967, 1.44825967, 1.44825967, 1.44825967, 1.44825967, 1.44825967, 1.448696, 1.448696, 1.448696, 1.448696, 1.44927778, 1.44852874, 1.44852874, 1.44852874, 1.44852874, 1.44852874, 1.44852874, 1.44852874, 1.44852874]

FRB_dec_rad = [0.57854271, 0.57830515, 0.57830515, 0.57830515, 0.57830515, 0.57830515, 0.57830515, 0.57834878, 0.57834878, 0.57834878, 0.57834878, 0.57834878, 0.57834878, 0.57825182, 0.57825182, 0.57826637, 0.57826637, 0.57747612, 0.57854029, 0.57854029, 0.57854029, 0.57854029, 0.57854029, 0.57854029, 0.57854029, 0.57854029]

FRB_time_window = np.ones_like(FRB_ra_rad)*args.dt/86400.

if args.dt > 1e3:
    n_bg_trials = args.bg_trials/1e2
    
    if args.dt > 1e5:
        n_sigs = np.r_[50:60:1, 60:100.1:2]
    elif args.dt <= 1e5:
        n_sigs = np.r_[2:10:1, 10:30.1:2]
        
elif args.dt <= 1e3:
    n_bg_trials = args.bg_trials
    n_sigs = np.r_[2:10:1, 10:30.1:2]
    
n_sig_trials = args.sig_trials

print("Starting Trials")

##### Defining functions to compute background trials, signal trials and flux ######    

def do_background_trials(N=n_bg_trials):
    src = cy.sources(FRB_ra_rad, FRB_dec_rad, mjd = FRB_mjd_time, sigma_t = np.zeros_like(FRB_ra_rad), t_100 = FRB_time_window)
    conf = {'extended':True, 'space':"ps",'time':"transient",'sig':"transient",'flux': cy.hyp.PowerLawFlux(args.gamma)}
    tr = cy.get_trial_runner(conf, src = src, ana=ana)
    trials = tr.get_many_fits(N,logging=False)
    dir = cy.utils.ensure_dir('{}/dt/{}'.format(bg_dir, args.dt))
    filename = '{}/bg_trials_seed_{}.npy'.format(dir, args.seed)
    print('->', filename)
    np.save(filename, trials.as_array)

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