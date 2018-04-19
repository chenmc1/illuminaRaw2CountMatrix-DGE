#!/bin/bash
##############################  DragonsTooth  dragonstooth1.arc.vt.edu
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=24
#PBS -W group_list=dragonstooth
#PBS -q normal_q
#PBS -A tomomics 
#PBS -e /work/dragonstooth/chenmc1/07_expTomics_2172/err.txt
#PBS -o /work/dragonstooth/chenmc1/07_expTomics_2172/out.txt
#PBS -M chenming.cui@icloud.com
#PBS -m bea

workingDir=/work/dragonstooth/chenmc1/07_expTomics_2172

###### this is the pipline final version. 03/29/2017

## this script used as pipeline
##################################
## 03232017 each accession take 20h~
## perspectus: this is dGE_analysis pipeline, designed to run in univac environment

## directory: 00_rawReads(start from unziped raw files), need manually make 00_directory and the zip files; 01_fastQC; 02_trimmedReads

##################################

#### load all the modules that necessary for the pipline from cluster
module load jdk/1.8.0
module load gcc/4.7.2
module load bowtie2
 

#### export all the PATHs for softwares
PATH=$HOME/blueridge/software/FastQC:$PATH
PATH=$HOME/blueridge/software/Trimmomatic-0.33:$PATH
PATH=$HOME/newriver/software/samtools-1.2:$HOME/newriver/software/star/bin/Linux_x86_64_static:$PATH
PATH=$HOME/newriver/software/trinityrnaseq-Trinity-v2.4.0/trinity-plugins/rsem-1.2.19/sam:$HOME/newriver/software/trinityrnaseq-Trinity-v2.4.0:$PATH
PATH=$PATH:/home/chenmc1/newriver/software/subread-1.5.0-p2-Linux-x86_64/bin



###
###
#####################  Step1. unzip rawreads.gz

for file in $workingDir/00_rawReads/*


    do

        if [[ $file =~ tgz$ ]]
        then

                tar -zxvf $file
                echo "file unziped"

        elif [[ $file =~ gz$ ]]
        then
                gunzip $file;
                echo "file unziped"
        else
                echo "there are files not in tgz or gz format that unable to unzip!"
        fi;

    done;
    



#:<<'COMMENT'
##################### Step2. quality control by fastqc

[ -d $workingDir/01_fastQC ] || mkdir $workingDir/01_fastQC


#for file in $workingDir/00_rawReads/*

 #   do

fastqc -f fastq -t 4 -o $workingDir/01_fastQC $workingDir/00_rawReads/*R1*.fastq $workingDir/00_rawReads/*R2*.fastq
   #    echo "FastQC is done successfully!"
    #done
    
COMMENT


#:<<'COMMENT'
###################### step3. triming #############

[ -d $workingDir/02_trimmedReads ] || mkdir $workingDir/02_trimmedReads
#mkdir $workingDir/02_trimmedReads
inputTrim=$workingDir/00_rawReads
outputTrim=$workingDir/02_trimmedReads
 
for file in $workingDir/00_rawReads/*_R1_001.fastq                
do                                                       
    withpath="${file}"
    filename="${withpath##*/}"  
    base="${filename%*_R1_001.fastq}"
    echo "${base}"                  
 
    java -jar $HOME/blueridge/software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 \
    $inputTrim/"${base}"_R1_001.fastq \
    $inputTrim/"${base}"_R2_001.fastq \
    $outputTrim/"${base}".R1.trim.fastq \
    $outputTrim/"${base}".R1.trimunpaired.fastq \
    $outputTrim/"${base}".R2.trim.fastq \
    $outputTrim/"${base}".R2.trimunpaired.fastq \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
done;

[ -d $outputTrim/trimUnpariedRaw ] || mkdir $outputTrim/trimUnpariedRaw

mv $outputTrim/*trimunpaired.fastq $outputTrim/trimUnpariedRaw

COMMENT




###
###
################################   Step4. Star alignment  #######################
###??? to index the genome,done! including the annotation for junction,,, just need the sj.out.tab file.
 
[ -d $workingDir/03_star ] || mkdir $workingDir/03_star

cd $workingDir/03_star
inputStar=$workingDir/02_trimmedReads
outputStar=$workingDir/03_star 

for file in $workingDir/02_trimmedReads/*.R1.trim.fastq                
                   
do                                                        
	withpath="${file}"
	filename="${withpath##*/}"    
	base="${filename%*.R1.trim.fastq}"  
	echo "${base}"  
	STAR --outReadsUnmapped Fastx --runThreadN 24 --genomeDir /home/chenmc1/newriver/reference_genome/tomato4star_3.0 \
	--readFilesIn $inputStar/"${base}".R1.trim.fastq $inputStar/"${base}".R2.trim.fastq --outFileNamePrefix $outputStar/"${base}".
 
 
    samtools view -Sb $outputStar/"${base}".Aligned.out.sam > $outputStar/"${base}".Aligned.out.bam
    # here actually no need to sort the bam file for space saving
    
    samtools sort $outputStar/"${base}".Aligned.out.bam $outputStar/"${base}".Aligned.sort.bam
done;

## next time :use outP.bam  outP.sort for the above



##############################
################################## here is star alignment for the trimmed unpaired end reads. use single end stratergy....

[ -d $workingDir/031_star_unpair ] || mkdir $workingDir/031_star_unpair
#mkdir $workingDir/031_star_unpair 

cd $workingDir/031_star_unpair
inputStarUn=$workingDir/02_trimmedReads/trimUnpariedRaw
outputStarUn=$workingDir/031_star_unpair 

for file in $workingDir/02_trimmedReads/*.R1.trim.fastq                
                   
do                                                        
	withpath="${file}"
	filename="${withpath##*/}"    
	base="${filename%*.R1.trim.fastq}"  
	echo "${base}"  
	STAR --outReadsUnmapped Fastx --runThreadN 24 --genomeDir /home/chenmc1/newriver/reference_genome/tomato4star \
	--readFilesIn $inputStarUn/"${base}".R1.trimunpaired.fastq --outFileNamePrefix $outputStarUn/"${base}".
 
 
    samtools view -Sb $outputStarUn/"${base}".Aligned.out.sam > $outputStarUn/"${base}".Aligned.out.bam
    
    samtools sort $outputStarUn/"${base}".Aligned.out.bam $outputStarUn/"${base}".Aligned.sort.bam
done;
####
##### use " outU.sort" outU.bam


:<<'COMMENT'
###############################  denovo assembly for the #03_star_unmapped_mate1/mate2 ####################

##add .fastq to each of the unmapped files

[ -d $workingDir/04_deno_trinity ] || mkdir $workingDir/04_deno_trinity
#mkdir $workingDir/04_deno_trinity
inputTri_de=$workingDir/03_star
outputTri_de=$workingDir/04_deno_trinity 


#for matefile in $inputTri_de/*Unmapped*
#do
#	mv "$matefile" "$matefile.fastq"
#done;

for file in $workingDir/02_trimmedReads/*.R1.trim.fastq                
                   
do                                                        
	withpath="${file}"
	filename="${withpath##*/}"    
	base="${filename%*.R1.trim.fastq}"  
	echo "${base}" 

    #mkdir $outputTri_de/"${base}"_trinity  
 
## add prefix for 
    Trinity --seqType fq --max_memory 50G --left $inputTri_de/"${base}".Unmapped.out.mate1.fastq --right $inputTri_de/"${base}".Unmapped.out.mate2.fastq \
--output $outputTri_de/"${base}"_trinity  
done;

######
######
######  ??
## this rename part can't run together with the previous one, need totally finish running the previous code





for file in $workingDir/02_trimmedReads/*.R1.trim.fastq                
                   
do                                                        
	withpath="${file}"
	filename="${withpath##*/}"    
	base="${filename%*.R1.trim.fastq}"  
	echo "${base}"
    mv $outputTri_de/"${base}"_trinity/Trinity.fasta $outputTri_de/"${base}"_trinity/"${base}"_Trinity.fasta
done;



####### merge bam #########


### test to merge the two bam file, didn't work now
workingDir=/work/newriver/chenmc1/06_DEG4000/01_expTomics_1589
for file in $workingDir/02_trimmedReads/*.R1.trim.fastq                
                   
do                                                        
	withpath="${file}"
	filename="${withpath##*/}"    
	base="${filename%*.R1.trim.fastq}"  
	echo "${base}"
	 
	samtools merge /work/newriver/chenmc1/06_DEG4000/01_expTomics_1589/03_star/"${base}".Aligned.outT.bam /work/newriver/chenmc1/06_DEG4000/01_expTomics_1589/03_star/"${base}".Aligned.out.bam /work/newriver/chenmc1/06_DEG4000/01_expTomics_1589/031_star_unpair/"${base}".Aligned.out.bam

done
COMMENT
###	




######### count reads###### feature count

gtfDir=/home/chenmc1/newriver/reference_genome/tomato4star_3.0
bamIn=$workingDir/03_star
countOut=$workingDir/05_readsCount 

[ -d $workingDir/05_readsCount ] || mkdir $workingDir/05_readsCount
#mkdir $workingDir/05_readsCount

for file in $workingDir/02_trimmedReads/*.R1.trim.fastq                
                   
do                                                        
	withpath="${file}"
	filename="${withpath##*/}"    
	base="${filename%*.R1.trim.fastq}"  
	echo "${base}"

    featureCounts -p -t exon -g gene_id -a  $gtfDir/ITAG3.10_gene_models_convertFromgff.gtf -o $countOut/"${base}".txt $bamIn/"${base}".Aligned.sort.bam
done

cat $workingDir/05_readsCount/*.txt > $workingDir/05_readsCount/rawcount.csv
# rename this total csv file and then copy the column infor file to a local dir.
# refer a python pipeline: https://github.com/maxplanck-ie/rna-seq-qc/blob/master/rna-seq-qc.py
