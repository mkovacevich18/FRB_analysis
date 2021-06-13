#!/bin/sh

base="/scratch/mkovacevich"
python_script_path="/data/user/mkovacevich/FRB_analysis"

rm -r "$base"/outputs
rm -r "$base"/errors
rm -r "$base"/job_files
rm -r "$base"/logs

mkdir "$base"/outputs/
mkdir "$base"/errors/
mkdir "$base"/logs/
mkdir "$base"/job_files/

mkdir "$base"/job_files/execs/

mkdir "$base"/job_files/subs/

mkdir "$base"/job_files/dags/

gammas=(2.0 3.0)
time_window=(0.01 0.1 1.0 10.0 100.0 1000.0 10000.0 100000.0 1000000.0 1000000.0)
choose=(1.0 2.0)

dag_path="$base/job_files/dags/FRB_121102_Catalog_trials.dag"
touch ${dag_path}

COUNTER=({1..100..1})
    
for gamma in ${gammas[@]}; do
    for dt in ${time_window[@]}; do
        for c in ${choose[@]}; do
	       for s in ${COUNTER[@]}; do 
                    exec_path="$base/job_files/execs/Job_gamma_${gamma}_dt_${dt}_seed_${s}_c${c}.sh"
		    touch ${exec_path}
		    echo "#!/bin/sh" >> ${exec_path} 
		    echo "python3 $python_script_path/FRB_121102_stacking_sens.py --gamma ${gamma} --dt ${dt} --seed ${s} --choose ${c} --bg_trials 10000.0 --sig_trials 100.0 --cpus 1.0" >> ${exec_path}
		    sub_path="$base/job_files/subs/Submit_gamma_${gamma}_dt_${dt}_seed${s}_c${c}.submit"
                    touch ${sub_path}
                    echo "executable = ${exec_path}" >> ${sub_path}
                    echo "output = $base/outputs/gamma_${gamma}_dt_${dt}_seed_${s}_c${c}.out" >> ${sub_path}
                    echo "error = $base/errors/gamma_${gamma}_dt_${dt}_seed_${s}_c${c}.err" >> ${sub_path}
                    echo "log = $base/logs/gamma_${gamma}_dt_${dt}_seed_${s}_c${c}.log" >> ${sub_path}
                    echo "getenv = true" >> ${sub_path}
                    echo "universe = vanilla" >> ${sub_path}
                    echo "notifications = never" >> ${sub_path}
                    echo "should_transfer_files = YES" >> ${sub_path}
                    echo "request_memory = 6GB" >> ${sub_path}
	            echo "request_cpus = 1" >> ${sub_path}
                    echo "queue 1" >> ${sub_path}
			
		    echo "JOB ${gamma}.${dt}.${s}.${c} ${sub_path}" >> ${dag_path}
		done

	done
   done
done
