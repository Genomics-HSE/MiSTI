#!/bin/bash


SAMPLE1=$1
SAMPLE2=$2
nthr=$3


# Функция для обработки одного региона
process_region() {
    local chr=$1
    local start=$2
    local end=$3
    local region="${chr}:${start}-${end}"
    local output_file="./2dsfs/${SAMPLE1}.${SAMPLE2}.${chr}_${start}_${end}.2dsfs"
    
    echo "Обрабатываю: $region"
    
    realSFS \
        "./${SAMPLE1}/${SAMPLE1}.sfs/${SAMPLE1}.all.saf.idx" \
        "./${SAMPLE2}/${SAMPLE2}.sfs/${SAMPLE2}.all.saf.idx" \
        -r "$region" > "$output_file"
}

# Экспорт функции и переменных
export -f process_region
export SAMPLE1 SAMPLE2

# Создание директории для результатов
mkdir -p ./2dsfs

# Параллельная обработка всех строк
cat window.grch38.bed | parallel --colsep '\t' --jobs ${nthr} --progress --joblog parallel.${SAMPLE1}.${SAMPLE2}.log \
    process_region {1} {2} {3}
