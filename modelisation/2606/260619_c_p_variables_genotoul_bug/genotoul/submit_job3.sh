#!/bin/bash
JOB2=38236497

while true; do
    COUNT=$(squeue -u $USER -h | wc -l)
    echo "$(date) - Jobs en queue : $COUNT"
    if [ $COUNT -le 1874 ]; then
        JOB3=$(sbatch --array=1-625 --dependency=afterany:$JOB2 \
               --export=ALL,TASK_OFFSET=2500 submit_c_p.sh | awk '{print $NF}')
        if [ -n "$JOB3" ]; then
            echo "Job 3 soumis : $JOB3"
            break
        else
            echo "Echec, on reessaie dans 5 min..."
        fi
    fi
    sleep 300
done
