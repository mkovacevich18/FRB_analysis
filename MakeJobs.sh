#!/bin/sh

#Executable file with code adapted from M. Campana's exectuable that creates jobs and dag files

mkdir '/scratch/mkovacevich/outputs/'
mkdir '/scratch/mkovacevich/errors/'
mkdir '/scratch/mkovacevich/logs/'
mkdir '/scratch/mkovacevich/job_files/'


mkdir '/scratch/mkovacevich/job_files/execs/'

mkdir '/scratch/mkovacevich/job_files/subs/'

mkdir '/scratch/mkovacevich/job_files/dags/'

gammas=(2.0 3.0)
time_window=(0.01 0.1 1.0 10.0 100.0 1000.0 10000.0 100000.0)
#Choose is equal to 1.0, 2.0 or 3.0 depending on whether bg trials (1.0), sig trials (2.0) or flux calculations (3.0) are being done
choose=(2.0)

dag_path="/scratch/mkovacevich/job_files/dags/Dagman_sig_trials.dag"
touch ${dag_path}

#Counter can change in order to optimize how many jobs can be run on npx, 1e9 bg trials were run for time windows less than or equal to 1000 seconds 1e7 bg trials
#1e7 bg trials are run for time windows greater than 1000 seconds
#The product of Counter * --bg_trials arguments below should produce the total number of bg trials needed for a time window

#10,000 signal trials were run for each time window, the produce of Counter * --sig_trials should produce 10,000
COUNTER=({1..50..1})
    
for gamma in ${gammas[@]}; do
    for dt in ${time_window[@]}; do
        for c in ${choose[@]}; do
	    #if [ c -ne 2.0 ]; then
	       for s in ${COUNTER[@]}; do 
                    exec_path="/scratch/mkovacevich/job_files/execs/Job_gamma_${gamma}_dt_${dt}_seed_${s}_c${c}.sh"
		    touch ${exec_path}
		    echo "#!/bin/sh" >> ${exec_path} 
		    echo "python3 /data/user/mkovacevich/FRB_analysis/FRB_stacking_sensitivity.py --gamma ${gamma} --dt ${dt} --seed ${s} --choose ${c} --bg_trials 20000.0 --sig_trials 100.0 --cpus 2.0" >> ${exec_path}
		    sub_path="/scratch/mkovacevich/job_files/subs/Submit_gamma_${gamma}_dt_${dt}_seed${s}_c${c}.submit"
                    touch ${sub_path}
                    echo "executable = ${exec_path}" >> ${sub_path}
                    echo "output = /scratch/mkovacevich/outputs/gamma_${gamma}_dt_${dt}_seed_${s}_c${c}.out" >> ${sub_path}
                    echo "error = /scratch/mkovacevich/errors/gamma_${gamma}_dt_${dt}_seed_${s}_c${c}.err" >> ${sub_path}
                    echo "log = /scratch/mkovacevich/logs/gamma_${gamma}_dt_${dt}_seed_${s}_c${c}.log" >> ${sub_path}
                    echo "getenv = true" >> ${sub_path}
                    echo "universe = vanilla" >> ${sub_path}
                    echo "notifications = never" >> ${sub_path}
                    echo "should_transfer_files = YES" >> ${sub_path}
                    echo "request_memory = 8GB" >> ${sub_path}
	            echo "request_cpus = 2" >> ${sub_path}
                    echo "queue 1" >> ${sub_path}
			
		    echo "JOB ${gamma}.${dt}.${s}.${c} ${sub_path}" >> ${dag_path}
		done
           #fi

	done
 done
done
