#!/bin/bash

#Mapping and MACS2 analysis ran for neo, old and CD4 T cells
#Daniel J. Simpson

module add trim_galore
module add bowtie2/2.5.0
module add samtools
module add picard/3.0.0
module add macs2
module load ucsc_tools/450


f="Inp_N1P2
Inp_N2P2
Inp_O1P2
Inp_O2P2
Inp_O3P2
IP_N1P2
IP_N2P2
IP_O1P2
IP_O2P2
IP_O3P2
EZH2_CD4_d32
EZH2_CD4_d43
EZH2_CD4_d48"

#Setting reference genome
bwt2_idx=/labs/vsebast/shared/refs/bw2/hg38/GRCh38_noalt_as

trim_galore --paired --gzip --cores 3 ${id}_*1.fq.gz ${id}_*2.fq.gz --output_dir ../Trimmed/

mkdir mapping_filtering

mkdir macs2_processing

cd mapping_filtering

echo "Start Mapping $f"

bowtie2 -X2000 -x $bwt2_idx -1 ../Trimmed/${f}_1_val_1.fq.gz -2 ../Trimmed/${f}_2_val_2.fq.gz 2> $f.log -p 18 -S ./${f}.sam

echo "End Mapping $f"

samtools view -F 1804 -f 2 -q 30 -u $f.sam -o $f.filtered -O SAM -@18
samtools sort -n $f.filtered -o $f.sorted -O SAM -@18
samtools fixmate -r $f.sorted $f.fixmate -O SAM -@18
samtools view -F 1804 -f 2 -u $f.fixmate -o $f.fixmate.filtered -O SAM -@18
samtools sort $f.fixmate.filtered -o $f.filtered.sorted -O SAM -@18

picard MarkDuplicates -I $f.filtered.sorted -O $f.marked.sam -METRICS_FILE $f.qc -VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --REMOVE_DUPLICATES false

samtools view -F 1804 -f 2 $f.marked.sam -o $f.dedup -O SAM -@18

rm $f.filtered $f.sorted $f.fixmate $f.fixmate.filtered $f.filtered.sorted $f.marked.sam




#running MACS on NEOs, OLDs and CD4s separately

cd ../macs2_processing/

Ns="N1P2 N2P2"

Os="O1P2 O2P2 O3P2"

CD4s="EZH2_CD4_d32 EZH2_CD4_d43 EZH2_CD4_d48"

for f in $Os; do

    macs2 predictd -i ../mapping_filtering/IP_$f.dedup

done


#Ns frags are 298 and 297 in that order. Using 297

#Os frags are 300, 299, and 298

#CD4 frags are 296 295 295


#NEO
FRAG=297
macs2 callpeak -t ../mapping_filtering/IP_N1P2.dedup ../mapping_filtering/IP_N2P2.dedup -c ../mapping_filtering/Inp_N1P2.dedup ../mapping_filtering/Inp_N2P2.dedup \
    -g hs -p 1e-2 --keep-dup all -B --nomodel --extsize $FRAG --SPMR -n NEO_P2_Merged

#Old
FRAG=299
macs2 callpeak -t ../mapping_filtering/IP_O1P2.dedup ../mapping_filtering/IP_O2P2.dedup ../mapping_filtering/IP_O3P2.dedup \
    -c ../mapping_filtering/Inp_O1P2.dedup ../mapping_filtering/Inp_O2P2.dedup ../mapping_filtering/Inp_O3P2.dedup\
    -g hs -p 1e-2 --keep-dup all -B --nomodel --extsize $FRAG --SPMR -n OLD_P2_Merged

#CD4
FRAG=295    
macs2 callpeak -t ../mapping_filtering/EZH2_CD4*.dedup \
        -c ../ENCODE_control/ENCFF918OQM.bam \
        -B --nomodel --SPMR -g hs -p 1e-2 --keep-dup all --extsize $FRAG \
        -n EZH2_CD4_pooled_ENCinp \
    ;



#leave uncommented which version is being ran
fileID="NEO_P2_Merged"
#fileID="OLD_P2_Merged"
#fileID="EZH2_CD4_pooled_ENCinp"

# FE first
macs2 bdgcmp -t ${fileID}_treat_pileup.bdg -c ${fileID}_control_lambda.bdg -o ${fileID}_fe.bdg -m FE  -p 0.00001

LC_COLLATE=C
sort -k1,1 -k2,2n ${fileID}_fe.bdg > ${fileID}_s_fe.bdg
awk '$1 !~ /_/ && $1 !~ /M/' ${fileID}_s_fe.bdg > ${fileID}_sf_fe.bdg
bedGraphToBigWig ${fileID}_sf_fe.bdg /labs/vsebast/DJS/genomes/HomoSap/UCSC_ver/hg38.chrom.sizes ${fileID}_fe.bw



####### Running peakcalling on individual files to get individual .narropeak files for DiffBind analysis later #############

list="N1P2
N2P2
O1P2
O2P2
O3P2"


FRAG=298


for fileID in $list; do

    macs2 callpeak -t ../mapping_filtering/IP_$fileID.dedup \
            -c ../mapping_filtering/Inp_$fileID.dedup \
            -g hs -p 1e-2 --keep-dup all -B --nomodel --extsize $FRAG --SPMR \
            -n ${fileID} \
            --outdir . \
        ;


done


#Running for CD4s

CD4s="EZH2_CD4_d32 
EZH2_CD4_d43 
EZH2_CD4_d48"

FRAG="295"

for fileID in $CD4s; do

  macs2 callpeak -t ../mapping_filtering/${fileID}.dedup \
          -c ../ENCODE_control/ENCFF918OQM.bam \
          -B --nomodel --SPMR -g hs -p 1e-2 --keep-dup all --extsize $FRAG \
          -n ${fileID} \
          --outdir . \
      ;


  macs2 bdgcmp -t ${fileID}_treat_pileup.bdg -c ${fileID}_control_lambda.bdg -o ${fileID}_fe.bdg -m FE  -p 0.00001


done