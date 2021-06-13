'''
Programmer: Mike Kovacevich
Last edited: 11/1/20
E-mail: mgk56@drexel.edu

Stacking sensitivity for different time windows and different gamma (spectral indices). Gamma will range from 2-3 and the time windows will range from 10^(-2) seconds to 10^5 seconds. Stacking sensitivity will be performed with Csky likelihood software.  
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import argparse
import histlite as hl
import csky as cy

#Building/loading MESC data from analysis directory
ana_dir = cy.utils.ensure_dir('/data/user/mkovacevich/FRB_analysis/cascades_ana')
repo = cy.selections.Repository()
ana = cy.analysis.Analysis(repo, cy.selections.MESEDataSpecs.mesc_7yr, dir=ana_dir)

#Directories to store the trials
trials_dir = cy.utils.ensure_dir('/data/user/mkovacevich/FRB_analysis/trials')
sig_dir = cy.utils.ensure_dir('{}/One_thousand_Sources_sig'.format(trials_dir))
bg_dir = cy.utils.ensure_dir('{}/One_thousand_Sources_bg'.format(trials_dir))

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

#Calculating sensitivities, 3sigma at 90% and discovery potential as functions of different time windows for 1 FRB
#The time windows will range from ~1 ms to 10^4 seconds (following previous times windows)
cy.CONF['ana'] = ana
cy.CONF['mp_cpus'] = args.cpus

#loading analysis object to get mjd of time data-set and events
a = ana.anas[0]

#good_indices represents an array that tracks the indices of FRB events that fall within livetime of MESC 7 yr
FRB_mjd_time = np.load('/data/user/mkovacevich/FRB_analysis/One_thousand_FRBS_MJD.npy', allow_pickle = True )
FRB_ra_rad = np.load('/data/user/mkovacevich/FRB_analysis/One_thousand_FRBS_RA.npy', allow_pickle = True)
FRB_dec_rad = np.load('/data/user/mkovacevich/FRB_analysis/One_thousand_FRBS_DEC.npy', allow_pickle = True)

FRB_time_window = np.ones_like(FRB_ra_rad)*args.dt/86400.

'''
#Now, we are going to simplify the above arrays to only include FRBs that burst during the livetime of MESC 7 yr
for i, time in enumerate(mjd_source_time):
    if float(time) >= min(a.data['mjd']) and float(time) <= max(a.data['mjd']):
        FRB_mjd_time.append(time)
        good_indices.append(i)
'''     

#Defining the number of signal and background trials to run
n_bg_trials = args.bg_trials
n_sig_trials = args.sig_trials

print("Starting Trials")

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
    # notice: trials.as_array is a numpy structured array, not a cy.utils.Arrays
    np.save(filename, trials.as_array)
    
n_sigs = np.r_[2:10:1, 10:30.1:2]
def do_signal_trials(n_sig, N=n_sig_trials):
    # get trial runner
    src = cy.sources(FRB_ra_rad, FRB_dec_rad, mjd = FRB_mjd_time, sigma_t = np.zeros_like(FRB_ra_rad), t_100 = FRB_time_window)
    conf = {'extended':True, 'space':"ps",'time':"transient",'sig':"transient",'flux': cy.hyp.PowerLawFlux(args.gamma)} 
    tr = cy.get_trial_runner(conf, src = src, ana=ana)
    # run trials
    trials = tr.get_many_fits(N, n_sig, poisson = True, logging=False)
    # save to disk
    dir = cy.utils.ensure_dir('{}/gamma/{}/dt/{}/n_sig/{}'.format(sig_dir, args.gamma, args.dt, n_sig))
    filename = '{}/sig_trials_{}.npy'.format(dir, args.seed)
    print('->', filename)
    # notice: trials.as_array is a numpy structured array, not a cy.utils.Arrays
    np.save(filename, trials.as_array)
    
def ndarray_to_TSD(trials):
    return cy.dists.TSD(cy.utils.Arrays(trials))

def find_n_sig(beta=0.9, nsigma=None):
    # get signal trials, background distribution, and trial runner
    sig_trials = cy.bk.get_best(sig, 'gamma', args.gamma, 'dt', args.dt, 'n_sig')
    b = cy.bk.get_best(bg, 'gamma', args.gamma, 'dt', args.dt)
    src = cy.sources(FRB_ra_rad, FRB_dec_rad, mjd = FRB_mjd_time, sigma_t = np.zeros_like(FRB_ra_rad), t_100 = FRB_time_window)
    conf = {'extended':True, 'space':"ps",'time':"transient",'sig':"transient",'flux': cy.hyp.PowerLawFlux(args.gamma)}
    tr = cy.get_trial_runner(conf, src = src, ana=ana)
    # determine ts threshold
    if nsigma is not None:
        ts = b.isf_nsigma(nsigma)
    else:
        ts = b.median()
    # include background trials in calculation
    trials = {0: b.trials}
    trials.update(sig_trials)
    # get number of signal events
    # (arguments prevent additional trials from being run)
    result = tr.find_n_sig(ts, beta, max_batch_size=0, logging=False, trials=trials, n_bootstrap=1)
    # return flux
    return tr.to_E2dNdE(result, E0=1e5)


#The bg and sig trials can be run at the same time used a diamond dag. Both of these trials need to be complete before
#going on to compute the sensitivities
if args.choose == 1:
    print('Computing bg trials w/ the time window ' + str(args.dt))
    do_background_trials(N=n_bg_trials)
    
elif args.choose == 2:
    print('Computing signal trials for gamma = ' + str(args.gamma) + 'and w/ the time window ' + str(args.dt))
    for n_sig in n_sigs:
        do_signal_trials(n_sig, N=n_sig_trials)
        
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