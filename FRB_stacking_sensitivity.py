'''
Programmer: Mike Kovacevich
Last edited: 11/22/20
E-mail: mgk56@drexel.edu

*****This script is used to only run the background and signal trials, FRB_stacking_sensitivity_flux_calculations.ipynb is used to calculate the fluxes*****

Stacking sensitivity for different time windows and different gamma (spectral indices). Gamma will range from 2-3 and the time windows 
will range from 10^(-2) seconds to 10^5 seconds and are log-spaced. Stacking sensitivity will be performed with Csky likelihood software.  
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
ana = cy.analysis.Analysis(repo, cy.selections.MESEDataSpecs.mesc_7yr, dir=ana_dir) # use load_sig = False if running background trials (reduces computational time)

#Directories to store the trials as nested dicts
trials_dir = cy.utils.ensure_dir('/data/user/mkovacevich/FRB_analysis/trials')
sig_dir = cy.utils.ensure_dir('{}/weighted_livetime_sig'.format(trials_dir))
bg_dir = cy.utils.ensure_dir('{}/weighted_livetime_bg'.format(trials_dir))

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

#Source times and equitorial coordinates of FRB single bursts
'''
source_times = ['2018-03-11T04:11:54.800','2018-03-09T02:49:32.990000','2018-03-01T07:34:19.760000','2017-12-09T20:34:23.500000','2017-09-22T11:22:23.400000','2017-08-27T16:20:18','2017-01-07T20:05:45.139000','2016-06-08T03:53:01.088000','2016-04-10T08:33:39.680000','2016-03-17T09:00:36.530000','2016-01-02T08:28:39.374000','2015-12-30T16:15:46.525000','2015-12-06T06:17:52.778000','2015-08-07T17:53:55.830000','2015-06-10T05:26:59.396000','2015-04-18T04:29:06.657000','2015-02-15T20:41:41.714000','2014-05-14T17:14:11.060000','2013-11-04T18:04:11.200000','2013-07-29T09:01:51.190000','2013-06-28T03:58:00.178000','2013-06-26T14:55:59.771000','2012-10-02T13:09:18.436000','2012-01-27T08:11:21.725000','2011-07-03T18:59:40.607000','2011-06-26T21:33:17.477000','2011-05-23T15:06:19.700000','2011-02-20T01:55:48.096000']
FRB_ra_deg = [322.89,321.18,93.18,237.60,322.46,12.33,170.69,114.17,130.35,118.45,339.70,145.21,290.35,340.63,161.11,109.15,274.36,338.52,101.04,205.34,135.76,246.77,273.70,348.77,352.71,315.93,326.30,338.66]
FRB_dec_deg = [-56.26,-32.02,4.56,-45.83,-6.01,-64.45,-4.98,-39.20,6.08,-28.39,-29.82,-2.55,-3.87,-54.92,-39.91,-18.99,-3.10,-11.69,-50.72,-4.0,3.44,-6.54,-84.80,-17.57,-1.13,-43.26,.016,-11.60]

#Converting date and times to mjd time
src_t = Time(source_times)
mjd_source_time = src_t.mjd
'''
#loading analysis object to get mjd of time data-set and events
a = ana.anas[0]

#good_indices represents an array that tracks the indices of FRB events that fall within livetime of MESC 7 yr
FRB_mjd_time = [57760.83732799769, 57547.161818148146, 57488.35670925926, 57464.37542280093, 57389.35323349537, 57386.67762181713, 57362.262416412035, 57241.74578506944, 57183.22707634259, 57130.18688260417, 57068.86228835648, 56791.718183564815, 56600.752907407405, 56502.376286921295, 56471.165279837965, 56469.62221957176, 56202.5481300463, 55953.341223668984, 55745.79144221065, 55738.89811894676, 55704.629394675925, 55612.080417777775]

FRB_ra_rad = [2.979102500229121, 1.992642407001926, 2.2750366799746087, 2.0673424989872835, 5.928883469024737, 2.534392606820966, 5.067563483165536, 5.945115031068284, 2.811899957888065, 1.9050268785518107, 4.788485335771643, 5.9082885838512045, 1.7634806762150708, 3.5838590860451562, 2.369458992507502, 4.306948995146407, 4.776966162708479, 6.087184832180623, 6.1559508047091995, 5.514018705825685, 5.695009349257497, 5.910732044803997]

FRB_dec_rad = [-0.08691739674931762, -0.6841690667817772, 0.10611601852125524, -0.4954989746411902, -0.5204571829447091, -0.0445058959258554, -0.06754424205218056, -0.9585348251952858, -0.6965609044709369, -0.33143802495372315, -0.054105206811824215, -0.2040289895581371, -0.8852309966115239, -0.06981317007977318, 0.060039326268604934, -0.11414453308042916, -1.4800392056911915, -0.3066543495754037, -0.01972222054753592, -0.755029434412747, 0.00027925268031909274, -0.20245819323134223]

FRB_time_window = np.ones_like(FRB_ra_rad)*args.dt/86400.

#Printing out length of each array to check they are all the same length, should each contain 22 elements
print(len(FRB_mjd_time))
print(len(FRB_ra_rad))
print(len(FRB_dec_rad))
print(len(FRB_time_window))

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

#This function is not typically run from this script since it does not need to run on npx, instead another script is run that utilizes this function
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


#The bg and sig trials can be run at the same time on npx. Both of these sig and bg trials need to be complete before calculating the fluxes
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

    
