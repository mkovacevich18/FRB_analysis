Analysis: Search for Coincident Neutrino Emission from Fast Radio Bursts (FRBs)

Location of scripts/notebooks (on cobalts): /data/user/mkovacevich/FRB_analysis/

Analysis Wiki: https://wiki.icecube.wisc.edu/index.php/A_Search_for_Coincident_Neutrino_Emission_from_Fast_Radio_Bursts_using_7_years_of_Cascades

Dataset: mese_cascades/version-001-p02

Analysis software: csky version 1.0.0 - the specific release of csky will be updated to the current main branch. Using the same seed, trials will be run to ensure there is no difference between my current version and the updated version

Before running scripts/notebooks on submitter/cobalt; the python version, meta-project and virtual environment are loaded (in that order). All bash scripts should contain the proper arguments in able to be run without any additional input.

Python version: /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/setup.sh
IceTray meta-project: /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-env combo/stable
Python virtual environment: source ~/py3venv/bin/activate

To run the bash scripts, the "base" and "python_script_path" file paths can be modified to write the job file/dagman/executable files to any directory


To generate signal and background trials, the following scripts/executables were used:
1. FRB_stacking_sensitivity.py
2. FRB_121102_stacking_sens.py

    1. These scripts have 7 arguments that are passed to them. 
        * gamma = spectral index used to inject events for signal trials
        * dt = duration of time window (10 time windows range from 0.01 s to 1e7 s and they are evenly log-spaced)
            - Most of this analysis only utilizes time windows up to 1e5 seconds, this analysis was extended past 1e5 on advice from the working group. Due to how many trials would need to be run, only the sensitivity/discovery potentials go out to 1e7 seconds. Bias tests, differential sensitivities and other time dependent quantities only go out to 1e5 seconds.
        * seed = this is can be used as a "seed" variable but is not in the trials, instead it tracks the total number of files that are made. For example, if we want to run 1e9 background trials, we can split this up into 100 or 1000 smaller jobs that can be run on submitter. Then the seed will go from [0,100] or [0, 1000]
        * choose = "choose" whether to run background trials (choose = 1.0), signal trials (choose = 2.0) or calculate sensitivity/discovery potentials (choose = 3.0). We should note that we only run the background and signal trials using this script. We calculate the sensitivity and discovery potentials using "FRB_stacking_sensitivity_flux_calculations.ipynb" since it is not computationally intensive to calculate these quantities and so we can use the cobalts. These scripts can generate the sensitivities/discovery potentials if necessary.
        * bg_trials = number of background trials to run in each job. For small time windows (dt < 1e3 seconds) 1e9 trials are used overall, for dt > 1e3 seconds we use 1e7 trials
            1. background trials are stored in '/data/user/mkovacevich/FRB_analysis/trials/bg'
        * sig_trials = number of signal trials to run in each job. Overall, we use 10,00 signal trials for all time windows
            1. signal trials are stored in '/data/user/mkovacevich/FRB_analysis/trials/sig'
        * cpus = number of cpus to use for each job
        
    2. For the first python script, the trials are stored in '/data/user/mkovacevich/FRB_analysis/trials/(bg or sig)'. For the second python script, the trials are stored in '/data/user/mkovacevich/FRB_analysis/FRB_121102_trials/(bg or sig)'
        
    3. To run on submitter using the following bash scripts: Sensitivity_MakeJobs.sh & FRB121102_Sensitivity_MakeJobs.sh
        1. These bash scripts are executable files that creates four directories that allow us to store information from the dagman file that is created. These executables create the dagman file but also creates the job files that run each python script with the arguments mentioned above.
        2. To run the dagman file, the following command is used "condor_submit_dag /scratch/mkovacevich/job_files/dags/(Insert Dagman File)" - only the dagman will be left after entering the dags directory so it is possible to tab complete the command. As mentioned earlier, the python version, metaproject and virtual environment are loaded before submitting jobs to submitter.
        
These two python scripts are identical except for the catalog that is used. The catalogs are hard-coded into each script; each catalog contains information about the source declination, right ascension and time (MJD). We assume an equal weighting for each source. The first script has a source list of 22 FRBs located in the southern sky; the second script is comprised of solely FRB121102 bursts. These scripts are run on submitter. This analysis intially focused on 1 catalog that corresponds to using the 1st script above. It was extended to FRB121101 on advice from the working group. Thus, most bias tests/other quantities were computed using the 1st catalog; this should not effect the analysis since only the source positions are effects, there is no additional weighting scheme.



To calculate the fluxes for sensitivities and discovery potentials, the following jupyter notebook was used:
1. FRB_stacking_sensitivity_flux_calculations.ipynb
    1. This python noteboook will compute the sensitivities and discovery potentials for both catalogs. These fluxes are then saved to /data/user/mkovacevich/FRB_analysis. The fluxes care computed using the pre-saved signal and background trials. This notebook is run on the cobalts since it is not *too* computationally intensive. The fluxes saved to the directory "Sensitivity_DiscoveryPotentials"
    2. This notebook is used to create the skymap plot that is located in the "FRB Sky Plots" section of my analysis wiki
    3. FRB121102_stacking_sensitivity.ipynb is used to get the RA, DEC, and MJD of each FRB121102 burst from the tns.csv file
    
    
    
To make plots of the fluxes, we use the following jupyter notebook:
1. Plots.ipynb
    1. This notebook will reproduce the plots shown in the "Sensitivities and Discovery Potentials" section of my analysis wiki. This notebook is run on the cobalts. From the saved fluxes, we calculate the energy-scaled flux per burst.
        * It should be noted that the saved fluxes are loaded relative to my FRB_analysis directory rather than the absolute file path
        
        
        
To test for any biases in the dataset, we use the following jupyter notebook:
1. BiasTest.ipynb
    1. This notebook creates plots that can be seen in "Bias Test" section of my analysis wiki. It computes the bias tests for gamma = 2, 3.
    2. Bias Tests were only performed on the first catalog. There is no major difference between the catalogs besides source position.
    
To plot the signal and background TS distributions, we use the following jupyter notebook:
1. Ninj_Pvalue.ipynb
    1. This notebook plots the overlaid background and signal TS distributions for gamma = 2 & 3 and for all time windows
    2. The plots in this notebook correspond to the plots in the section "Background and Signal TS Distributions" of my analysis wiki
    
    
    
To observe the angular resolution of the MESC dataset, we use the following jupyter notebook:
1. AnalysisProperties.ipynb
    1. This notebook serves to soley make a plot of delta_psi vs. energy. It is not in my analysis wiki.
    
    
    
To calculate the differential sensitivity, we use the following python script and jupyter notebook:
1. DifferentialSensitivity.py
    1. This script is run on submitter and computes the differential sensitivity. It then saves the output. The files have been moved to the directory "Differential_Sensitivity_trials".
    2. To run this script, we use the following bash script: "MakeJobs.sh"
    3. To make plots of the differential sensitivity, we use the notebook: "DifferentialSensitivity.ipynb"
        * This notebook will replicate the 4 plots that are shown in the sections: "Differential Sensitivity & Model Spectra" and "Differential Sensitivity Checks" of my analysis wiki
        1. The model spectra are digitized from "Neutrino Counterparts of FRBs" and stored in "Extracted_MMF_sensitivity.ipynb"
    4. To calculate the central 90% energy range, we use the script "Reduced_Energy_Range.py". This script will inject events over various energy ranges in order to observe which sensitivities correspond to the centeral 90% energy range for each time window.
        1. To actually calculate the sensitivity for each energy range, we use the notebook "Find_Reduced_Energy_Range.ipynb"
        2. To run this python script on submitter, the bash script "Reduced_Energy_Range.sh" is used (it contains all the appropriate arguments)
        3. This script only requires signal trials to be run, they are stored in '/data/user/mkovacevich/FRB_analysis/trials/reduced_energy_range_sig'
            
To produce the plot in the "Future Plans" we use the following jupyter notebook: 
1. N_Source_Stacking_Sensitivitiy.ipynb
    1. This notebook utilizes the sources that have been generated with the notebook "Generate_isotropic_FRBs.ipynb" in order to make a plot of how the sensitivity per burst decreases as more FRBs are observed. The generated sources are saved and loaded as "xxxx_FRBs_yyyy.npy" where xxxx = Ten, One_hundred, or One_thousand and yyyy = RA, DEC, or MJD
    2. To calculate the sensitivities, we need to run signal and background trials for the various number of sources. The signal and background trials are stored in: '/data/user/mkovacevich/FRB_analysis/trials/xxxx_Sources_bg', '/data/user/mkovacevich/FRB_analysis/trials/xxxx_Sources_sig' where xxxx = Ten, One_hundred, or One_thousand
    
To unblind the data, we will use the following python script and run on submitter:
1. Unblind.py
    1. This script can be run on cobalt. 
    2. To run the script, we must give it the following arguments that are listed below
        * catalog = 1.0 or 2.0 (1.0 -> main catalog, 2.0 -> FRB121102 catalog)
        * choose = 1.0 or 2.0 (1.0 -> calculate unblinded pre-trial pvalues, 2.0 -> calculate unblinded post-trial pvalues)
        Example of running script: "python3 Unblind.py --catalog 1.0 --choose 1.0"
    3. The purpose of this script is to obtain our unblinded pre and post trial p-values. After computing the pre-trial p-values, we store them as a dictionary that can be loaded when computing the post-trial results. To obtain the post-trial results, we will utilize correlated background trials and take the "hottest" p-value from each trial. This will form a distribution and we can see where the unblinded pre-trial p-values falls in this distribution. From this, we can obtain the post-trial p-value
            

        

    
    







