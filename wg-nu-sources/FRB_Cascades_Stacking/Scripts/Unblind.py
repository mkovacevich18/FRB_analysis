import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import argparse
import histlite as hl
import csky as cy

ana_dir = cy.utils.ensure_dir('/data/user/mkovacevich/FRB_analysis/cascades_ana')
repo = cy.selections.Repository()
ana = cy.analysis.Analysis(repo, cy.selections.MESEDataSpecs.mesc_7yr, dir=ana_dir)

parser = argparse.ArgumentParser(description='Choose the powerlaw flux and which catalog to unblind')
parser.add_argument('--catalog',type=float,help='choose which catalog to unblind')
parser.add_argument('--choose',type=float,help='choose between calculating pre and post trial results (pvalues)')
args = parser.parse_args()

cy.CONF['ana'] = ana
cy.CONF['mp_cpus'] = 10.0

bg_trials = {}

if args.catalog = 1.0:
    #MJD, RA and Dec of each FRB
    FRB_mjd_time = [57760.83732799769, 57547.161818148146, 57488.35670925926, 57464.37542280093, 57389.35323349537, 57386.67762181713, 57362.262416412035, 57241.74578506944, 57183.22707634259, 57130.18688260417, 57068.86228835648, 56791.718183564815, 56600.752907407405, 56502.376286921295, 56471.165279837965, 56469.62221957176, 56202.5481300463, 55953.341223668984, 55745.79144221065, 55738.89811894676, 55704.629394675925, 55612.080417777775]

    FRB_ra_rad = [2.979102500229121, 1.992642407001926, 2.2750366799746087, 2.0673424989872835, 5.928883469024737, 2.534392606820966, 5.067563483165536, 5.945115031068284, 2.811899957888065, 1.9050268785518107, 4.788485335771643, 5.9082885838512045, 1.7634806762150708, 3.5838590860451562, 2.369458992507502, 4.306948995146407, 4.776966162708479, 6.087184832180623, 6.1559508047091995, 5.514018705825685, 5.695009349257497, 5.910732044803997]

    FRB_dec_rad = [-0.08691739674931762, -0.6841690667817772, 0.10611601852125524, -0.4954989746411902, -0.5204571829447091, -0.0445058959258554, -0.06754424205218056, -0.9585348251952858, -0.6965609044709369, -0.33143802495372315, -0.054105206811824215, -0.2040289895581371, -0.8852309966115239, -0.06981317007977318, 0.060039326268604934, -0.11414453308042916, -1.4800392056911915, -0.3066543495754037, -0.01972222054753592, -0.755029434412747, 0.00027925268031909274, -0.20245819323134223]
    
#FRB121102 catalog    
elif args.catalog = 2.0:
    FRB_mjd_time = [57623.74402686, 57364.20463266, 57345.4624131, 57345.4575953, 57345.45248793, 57345.44769125, 57339.35604601, 57175.74828727, 57175.74762485, 57175.74566583, 57175.74351039, 57175.74283934, 57175.7425767, 57175.69972782, 57175.69314323, 57159.74422362, 57159.73760083, 56233.27492181, 57633.67986367, 57633.69515937, 57638.49937435, 57643.45730263, 57645.42958602, 57646.4660065, 57648.4369149 , 57649.45175697]

    FRB_ra_rad = [1.44852874, 1.44847784, 1.44847784, 1.44847784, 1.44847784, 1.44847784, 1.44847784, 1.44825967, 1.44825967, 1.44825967, 1.44825967, 1.44825967, 1.44825967, 1.448696, 1.448696, 1.448696, 1.448696, 1.44927778, 1.44852874, 1.44852874, 1.44852874, 1.44852874, 1.44852874, 1.44852874, 1.44852874, 1.44852874]

    FRB_dec_rad = [0.57854271, 0.57830515, 0.57830515, 0.57830515, 0.57830515, 0.57830515, 0.57830515, 0.57834878, 0.57834878, 0.57834878, 0.57834878, 0.57834878, 0.57834878, 0.57825182, 0.57825182, 0.57826637, 0.57826637, 0.57747612, 0.57854029, 0.57854029, 0.57854029, 0.57854029, 0.57854029, 0.57854029, 0.57854029, 0.57854029]

#Loading functions to load bg trials
def ndarray_to_TSD(trials):
    return cy.dists.TSD(cy.utils.Arrays(trials))

def tsd_merge(x):
    ts_values = np.concatenate([xx[0] for xx in x])
    n_zero = sum(xx[1] for xx in x)
    return cy.dists.TSD(ts_values, n_zero=n_zero)

#Length of the time window after a burst in which to search for neutrinos
FRB_time_window = np.ones_like(FRB_ra_rad)*args.dt/86400.

#loading smallest and largest time windows bg trials
time_windows = np.logspace(-2, 7, 10)

if args.chooose == 1.0:
    #Storing best fit values, ts, gamma and time window
    unblinding_results = {"time window" : [], "ts" : [], "ns" : [], "gamma" : [], "flux" : [], "pvalue" : [], "sigma"}


    #Loading bg trials to establish background expectations
    for dt in time_windows:
        print(dt)

        bg_trials[dt] = cy.bk.get_all(
        '/data/user/mkovacevich/FRB_analysis/trials/bg/dt/{}/'.format(dt),
        'bg_trials_seed_*.0.npy',
        pre_convert=lambda x: (x['ts'][x['ts'] > 0], np.sum(x['ts'] == 0)),
        merge = tsd_merge,
        log=True)
        bg = bg_trials[dt]

        FRB_time_window = np.ones_like(FRB_ra_rad)*dt/86400.

        src = cy.sources(FRB_ra_rad, FRB_dec_rad, mjd = FRB_mjd_time, sigma_t = np.zeros_like(FRB_ra_rad), t_100 = FRB_time_window)    
        conf = {'extended':True, 'space':"ps",'time':"transient",'sig':"transient"}
        tr = cy.get_trial_runner(conf, src = src, ana=ana)
        ts, ns, gamma = tr.get_one_fit(seed = 1, TRUTH=False)  #TRUTH = True for actual unblinding
        flux = tr.to_E2dNdE(ns, gamma=gamma E0=1e5)
        pvalue = bg.sf(ts)
        sigma = bg.sf_nsigma(ts)

        #pre-trial results, pvalue below needs the post-trial correction
        unblinding_results["time window"].append(dt)
        unblinding_results["ts"].append(ts)
        unblinding_results["ns"].append(ns)
        unblinding_results["gamma"].append(gamma)
        unblinding_results["flux"].append(flux)
        unblinding_results["pvalue"].append(pvalue)
        unblinding_results["sigma"].append(sigma)

    np.save('/data/user/mkovacevich/FRB_analysis/Unblinded_Results/unblinding_results.npy', unblinding_results, allow_pickle = True)
    
#Runnng correlated trials to get post trial results
elif args.choose == 2.0

    unblinding_results = np.load('/data/user/mkovacevich/FRB_analysis/Unblinded_Results/unblinding_results.npy', allow_pickle = True)
    index = np.linspace(0.0, 9.0, num=10)
    
    for (i, dt) in zip(index, time_windows):
        
        FRB_time_window = np.ones_like(FRB_ra_rad)*dt/86400.
        src = cy.sources(FRB_ra_rad, FRB_dec_rad, mjd = FRB_mjd_time, sigma_t = np.zeros_like(FRB_ra_rad), t_100 = FRB_time_window)
        conf = {'extended':True, 'space':"ps",'time':"transient",'sig':"transient"}
        
        #Accouting for correlated background trials (if pre-trial TS values are non-zero)
        trs = [cy.get_trial_runner(src=s) for s in src]
        tr_inj = cy.get_trial_runner(src=src, inj_conf=dict(src=src))

        #May need to increase the number of trials, bg for each source
        bgs = [cy.dists.TSD(tr.get_many_fits(1e3)) for tr in trs]

        #Setting up the correlated hypothesis
        multr = cy.trial.MultiTrialRunner(ana, tr_inj, trs, bgs = bgs,  mp_cpus=cy.CONF['mp_cpus'])

        trials = multr.get_many_fits(1000, n_sig=10, seed=1)

        bg_trials = multr.get_many_trials(10e3, seed=1)
        bg_trials['sum_ts'] = np.sum([bg_trials[name] for name in bg_trials.keys() if 'ts' in name], axis=0)
        
        #Getting the hottest p-value from each trial
        mlog10p = np.array([bg_trials[name] for name in bg_trials.keys() if name.startswith('mlog10p_')])
        ts = np.array([bg_trials[name] for name in bg_trials.keys() if name.startswith('ts_')])
        ns = np.array([bg_trials[name] for name in bg_trials.keys() if name.startswith('ns_')])
        gamma = np.array([bg_trials[name] for name in bg_trials.keys() if name.startswith('gamma_')])
        
        imax = np.argmax(mlog10p, axis=0)
        r = np.arange(len(imax))
        bg_trials_hottest = cy.utils.Arrays(dict(mlog10p=mlog10p[imax,r], ts=ts[imax,r], ns=ns[imax,r], gamma=gamma[imax,r]))
        
        #Load pre-trial values dictionary and compare the correlated bg trials to see where the unblinded ts/pvalues fall
        ts = unblinding_results["ts"][i]
        post_trial_pvalue = bg_trials_hottest.sf(ts)
        post_trial_sigma = bg_trials_hottest.sf_nsigma(ts)
        
        unblinding_results["post trial p-value"].append(post_trial_pvalue)
        unblinding_results["post_trial_sigma"].append(post_trial_sigma)
        
    np.save('/data/user/mkovacevich/FRB_analysis/Unblinded_Results/unblinding_results.npy', unblinding_results, allow_pickle = True)
        
        
        
    
    
    




    

    
