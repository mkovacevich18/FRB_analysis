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

dag_path="$base/job_files/dags/Differential_Sensitivity_Trials.dag"
touch ${dag_path}
    
for gamma in ${gammas[@]}; do
    for dt in ${time_window[@]}; do
                    exec_path="$base/job_files/execs/Job_gamma_${gamma}_dt_${dt}.sh"
		    touch ${exec_path}
		    echo "#!/bin/sh" >> ${exec_path} 
		    echo "python3 $python_script_path/DifferentialSensitivity.py --gamma ${gamma} --dt ${dt}  --cpus 2.0" >> ${exec_path}
		    sub_path="$base/job_files/subs/Submit_gamma_${gamma}_dt_${dt}.submit"
                    touch ${sub_path}
                    echo "executable = ${exec_path}" >> ${sub_path}
                    echo "output = $base/outputs/gamma_${gamma}_dt_${dt}.out" >> ${sub_path}
                    echo "error = $base/errors/gamma_${gamma}_dt_${dt}.err" >> ${sub_path}
                    echo "log = $base/logs/gamma_${gamma}_dt_${dt}.log" >> ${sub_path}
                    echo "getenv = true" >> ${sub_path}
                    echo "universe = vanilla" >> ${sub_path}
                    echo "notifications = never" >> ${sub_path}
                    echo "should_transfer_files = YES" >> ${sub_path}
                    echo "request_memory = 5GB" >> ${sub_path}
	            echo "request_cpus = 2" >> ${sub_path}
                    echo "queue 1" >> ${sub_path}
			
		    echo "JOB ${gamma}.${dt} ${sub_path}" >> ${dag_path}
   done
done
