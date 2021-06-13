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

gammas=(2.0)
time_window=(0.01 0.1 1.0 10.0 100.0 1000.0 10000.0 100000.0)
energy_exponent=(2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8)
choose=(2.0)

dag_path="$base/job_files/dags/Reduced_Energy_Range_Trials.dag"
touch ${dag_path}

COUNTER=({1..10..1})
    
for gamma in ${gammas[@]}; do
    for dt in ${time_window[@]}; do
        for c in ${choose[@]}; do
            for energy in ${energy_exponent[@]}; do
               for s in ${COUNTER[@]}; do 
                        exec_path="$base/job_files/execs/Job_gamma_${gamma}_dt_${dt}_seed_${s}_c${c}_energy${energy}.sh"
                touch ${exec_path}
                echo "#!/bin/sh" >> ${exec_path} 
                echo "python3 $python_script_path/Reduced_Energy_Range.py --gamma ${gamma} --dt ${dt} --seed ${s} --choose ${c} --bg_trials 10000.0 --sig_trials 1000.0 --cpus 1.0 --energy_exponent ${energy}" >> ${exec_path}
                sub_path="$base/job_files/subs/Submit_gamma_${gamma}_dt_${dt}_seed${s}_c${c}_energy${energy}.submit"
                        touch ${sub_path}
                        echo "executable = ${exec_path}" >> ${sub_path}
                        echo "output = $base/outputs/gamma_${gamma}_dt_${dt}_seed_${s}_c${c}_energy${energy}.out" >> ${sub_path}
                        echo "error = $base/errors/gamma_${gamma}_dt_${dt}_seed_${s}_c${c}_energy${energy}.err" >> ${sub_path}
                        echo "log = $base/logs/gamma_${gamma}_dt_${dt}_seed_${s}_c${c}_energy${energy}.log" >> ${sub_path}
                        echo "getenv = true" >> ${sub_path}
                        echo "universe = vanilla" >> ${sub_path}
                        echo "notifications = never" >> ${sub_path}
                        echo "should_transfer_files = YES" >> ${sub_path}
                        echo "request_memory = 5GB" >> ${sub_path}
                        echo "request_cpus = 1" >> ${sub_path}
                        echo "queue 1" >> ${sub_path}

                echo "JOB ${gamma}.${dt}.${s}.${c}.${energy} ${sub_path}" >> ${dag_path}
            done
        done
	done
   done
done
