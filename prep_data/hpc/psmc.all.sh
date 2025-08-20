#!/bin/bash



SAMPLE=$1
#SBATCH --job-name=psmc.sequential_jobs
#SBATCH --output=sequential_%j.log

#SBATCH --partition=compute

# Определяем количество ядер для каждого скрипта
N1=2
N2=24
N3=8
N4=1


# Массив с именами скриптов
SCRIPTS=("./coverage.sh ${SAMPLE} 2" "./split.cram.sh ${SAMPLE} 8" "./make.fq.sh ${SAMPLE} 8" "./psmc.sh ${SAMPLE}" )
CORES=($N1 $N2 $N3 $N4 )

# Запускаем первый скрипт
echo "Запуск первого скрипта с $N1 ядрами..."
JOB1=$(sbatch -A proj_1669 --parsable --cpus-per-task=$N1 ${SCRIPTS[0]})

# Последующие скрипты зависят от завершения предыдущих
JOB2=$(sbatch -A proj_1669 --parsable --dependency=afterok:$JOB1 --cpus-per-task=$N2 ${SCRIPTS[1]})
JOB3=$(sbatch -A proj_1669 --parsable --dependency=afterok:$JOB2 --cpus-per-task=$N3 ${SCRIPTS[2]})
JOB4=$(sbatch -A proj_1669 --parsable --dependency=afterok:$JOB3 --cpus-per-task=$N4 ${SCRIPTS[3]})


echo "Все задачи отправлены в очередь:"
echo "Job1 ID: $JOB1"
echo "Job2 ID: $JOB2" 
echo "Job3 ID: $JOB3"
echo "Job4 ID: $JOB4"

