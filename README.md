# Plant-LncRNA-pipeline-v2

![LncRNA pipline-V2](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2/blob/master/Plant-LncPipeV2.workflow.png) 





# **Plant-LncPipe-v2**





### **A pipeline for identifying and characterizing lncRNAs in plants.**

#### **Email:** xuechan.tian@bjfu.edu.cn;  jianfeng.mao@umu.se








## **1. Dependencies** 



    HISAT2
    StringTie
    LncFinder
    CPAT
    Diamond





## **2.Installation**



### **2.1. Install HISAT2**
    git clone https://github.com/DaehwanKimLab/hisat2.git
    cd hisat2
    make



### **2.2. Install StringTie**
    wget https://github.com/gpertea/stringtie/releases/download/v2.1.4/stringtie-2.1.4.Linux_x86_64.tar.gz
    tar xzf stringtie-2.1.4.Linux_x86_64.tar.gz
    export PATH=$PATH:/path/stringtie-2.1.6.Linux_x86_64/

	
### **2.3. Install PlantLncBoost**
    git clone https://github.com/xuechantian/PlantLncBoost.git


	
### **2.4. Install CPAT-plant and LncFinder-plant**
    git clone https://github.com/xuechantian/Plant-LncRNA-pipline.git


CPAT (version 1.2.4)


    conda create -n py27 python=2.7 -y
    source activate py27
    pip2 install CPAT


 
LncFinder-plant R Package

 
    install.packages("LncFinder")
    install.package("seqinr")
	
	
	
### **2.5. Install Diamond**
     conda install -c bioconda diamond
     
     

### **2.6. Install FEELnc**
    git clone https://github.com/tderrien/FEELnc.git
    export FEELNCPATH=/path/FEELnc/bin/
    export PERL5LIB=$PERL5LIB:/path/FEELnc/lib/
    export PATH=$PATH:/path/FEELnc/scripts/

	
	
### **2.7. Install fastp**
    conda install -c bioconda fastp	
	
	
	
## **3. Download Example Data (Soybean Data as an Example)**



### **3.1. Genome Assembly and Annotation File**
    Download the genome assembly fasta file and gff file from Phytozome 13 for Glycine max Wm82.a6.v1


### **3.2. Download Transcriptome Sequencing Data**
Example content of sra.list file


    SRR1174214
    SRR1174217
    SRR1174218
    SRR1174218
    SRR1174232
    ......


 
Download


    prefetch --option-file sra.list
	
	

 
## **4. Convert sra to fastq**
    fastq-dump SRR1174214.sra
    fastq-dump SRR1174217.sra
    fastq-dump SRR1174218.sra
    fastq-dump SRR1174228.sra
    fastq-dump SRR1174232.sra

	
	
## **5. Data Filtering and Quality Control**
    fastp -i SRR1174214.fastq -o SRR1174214_clean.fastq 
    fastp -i SRR1174217.fastq -o SRR1174217_clean.fastq 
    fastp -i SRR1174218.fastq -o SRR1174218_clean.fastq 
    fastp -i SRR1174228.fastq -o SRR1174228_clean.fastq 
    fastp -i SRR1174232.fastq -o SRR1174232_clean.fastq 


## **6. Run HISAT2 to map RNA-seq reads to the reference genome**



### **6.1. Construct reference genome**

    hisat2-build -p 8 genome.fasta genome.index 

	

### **6.2. Genome alignment with hisat2**
#### **Single-End RNA-seq  (e.g., Soybean)**


If the RNA-seq library is strand-specific, add the parameter "--rna-strandness RF".


    for i in `cat sra.list`; do hisat2 --new-summary --rna-strandness RF -p 10 -x genome.index ${i}_clean.fastq -S ${i}.sam; done

 
If the RNA-seq library is not strand-specific, remove  the parameter "--rna-strandness RF".


    for i in `cat sra.list`; do hisat2 --new-summary -p 10 -x genome.index ${i}_clean.fastq -S ${i}.sam; done

	
#### **Paired-End RNA-seq**


 If the RNA-seq library is strand-specific, add the parameter "--rna-strandness RF".

 
    for i in `cat sra.list`; do hisat2 --new-summary --rna-strandness RF -p 10 -x genome.index -1 ${i}_1_clean.fastq -2 ${i}_2_clean.fastq -S${i}.sam; done
	
 If the RNA-seq library is not strand-specific, remove  the parameter "--rna-strandness RF".

 
    for i in `cat sra.list`; do hisat2 --new-summary -p 10 -x genome.index -1 ${i}_1_clean.fastq -2 ${i}_2_clean.fastq -S ${i}.sam; done

	
	
### **6.3. Sort and compress sam files with samtools**

    for i in `cat sra.list`; do samtools view -S -b ${i}.sam | samtools sort -o ${i}.bam; done
	

	

## **7. Assemble transcripts using StringTie**

### **7.1. Format of Glycine_max_longest.gtf**
    Gm01	phytozomev13	exon 	103572	103594	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2" 
    Gm01	phytozomev13	CDS 	103572	103594	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2"
    Gm01	phytozomev13	exon 	103222	103288	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2"	
    Gm01	phytozomev13	CDS 	103222	103288	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2"	
    Gm01	phytozomev13	exon	102790	102857	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2"	
    Gm01	phytozomev13	CDS	102790	102857	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2";	
    Gm01	phytozomev13	exon	78986	79111	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2";	
    Gm01	phytozomev13	CDS	78986	79111	transcript_id "Glyma.01G000100.2"; gene_id "Glyma.01G000100.2";	

	
	
### **7.2. Transcription reconstruction of a single sample**

If strand-specific, add the parameter "--rf".


	for i in `cat sra.list`; do stringtie -p 10 --rf -G Glycine_max_longest.gtf  -o ${i}.gtf  ${i}.bam; done

 
If not strand-specific, remove the parameter "--rf".


    for i in `cat sra.list`; do stringtie  -p 10 -G Glycine_max_longest.gtf -o ${i}.gtf  ${i}.bam; done
	
	
### **7.3. Merge transcripts from multiple samples**	
	stringtie --merge -o merge.gtf  -G Glycine_max_longest.gtf  SRR*.gtf
	grep 'transcript_id "MSTRG' merge.gtf > candidate_transcript.gtf
	gffread -w candidate_transcript.fasta -g genome.fasta candidate_transcript.gtf
	grep '>' candidate_transcript.fasta | awk '{print \$1}' | sed 's/>//g' | sort -u > candidate_transcript.txt

	
	
## **8. LncRNA identification**	



### **8.1. Remove transcripts shorter than 200 bp and overlapping with known mRNAs**

    FEELnc_filter.pl -i candidate_transcript.gtf -a Glycine_max_longest.gtf --monoex=-1 -s 200 -p 20 > candidate_lncRNA.gtf
    cut -d ";" -f 2 candidate_lncRNA.gtf |sed 's/ transcript_id //g' | sed 's/"//g' | sort -u > candidate_lncRNA.txt


	
### **8.2. Identification of lncRNA with PlantLncBoost**	
#### **8.2.1. Dependencies**
    Python (>=3.7.3)
    Biopython
    NumPy
    Pandas
    SciPy
    CatBoost


#### **8.2.2. run PlantLncBoost**

#### **Feature extraction**
    python PlantLncBoost/Script/Feature_extraction.py -i candidate_transcript.fasta -o PlantLncBoost_feature.csv
	
#### **LncRNA prediction**
In the second column (Predicted_label) of the result file, 1 represents lncRNA and 0 represents mRNA.

    python PlantLncBoost/Script/PlantLncBoost_prediction.py -i PlantLncBoost_feature.csv -m ./PlantLncBoost/Model/PlantLncBoost_model.cb -t 0.5 -o PlantLncBoost_prediction.csv

	
### **8.3.  Identification of lncRNA with LncFinder-plant**	

R Package

    R
    library(LncFinder)
    library(seqinr)
	
	
import training data

    mRNA <- seqinr::read.fasta(file ="./data/training/mRNA.fasta")
    lncRNA <- seqinr::read.fasta(file ="./data/training/lncRNA.fasta")
    
    
Use "make_frequencies" function to generate the feature file

    frequencies <- make_frequencies(cds.seq = mRNA, lncRNA.seq = lncRNA, SS.features = FALSE, cds.format = "DNA", lnc.format = "DNA", check.cds = TRUE, ignore.illegal = TRUE)	
	
	
loading the model

    plant = readRDS("./Model/Plant_model.rda")
	
	
Identification of lncRNA 

    Seqs <- seqinr::read.fasta(file ="candidate_transcript.fasta")
    Plant_results <- LncFinder::lnc_finder(Seqs, SS.features = FALSE, format = "DNA", frequencies.file = frequencies, svm.model = plant, parallel.cores = 2)
	
	
Export results

    write.table(Plant_results, file ="plant-lncFinder.txt", sep ="\t",row.names =TRUE, col.names =TRUE,quote =FALSE)




### **8.4. Identification of lncRNA with CPAT-plant**
The coding probability (CP) cutoff: 0.46 (CP >=0.46 indicates coding sequence, CP < 0.46 indicates noncoding sequence).

    source activate py27
    cpat.py -x ./Model/Plant_Hexamer.tsv -d ./Model/Plant.logit.RData -g candidate_transcript.fasta -o CPAT_plant.output





### **8.5. Alignment of sequences to the UniProt protein database with diamond**
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gunzip uniprot_sprot.fasta.gz
    diamond makedb --in uniprot_sprot.fasta -d uniprot_out
    diamond blastx -d uniprot_out -q candidate_transcript.fasta -o uniprotoutput.txt




### **8.6. By intersecting the results obtained from the aforementioned steps, a set of high-confidence lncRNAs were obtained**

The id of the lncRNA

    Rscript prediction_insersection.sh candidate_lncRNA.txt PlantLncBoost_prediction.csv CPAT_plant.output plant-lncFinder.txt uniprotoutput.txt
    
The lncRNA gtf file

    grep -Fwf Final_lncRNA_results.txt candidate_transcript.gtf > lncRNA.gtf

	

  
	
## **9. Classify the final set of lncRNAs based on their genomic locations and sequence features**
Classification result file

	FEELnc_classifier.pl -i lncRNA.gtf -a Glycine_max_longest.gtf > lncRNA_classes.txt
	
	
Antisense_exonic-lncRNA

	awk -F '\t' '{if($1==1 && $6 == "antisense" && $10 == "exonic") {print $0}}' lncRNA_classes.txt > LncRNA_antisense_exonic.txt


Intronic-lncRNA

	awk -F '\t' '{if($1==1 && $6 == "sense" && $10 == "intronic") {print $0}}' lncRNA_classes.txt > LncRNA_intronic.txt


Upstream-lncRNA

	awk -F '\t' '{if($1==1 && $6 == "sense" && $7 == "intergenic" && $8 <= 2000 && $10 == "upstream") {print $0}}' lncRNA_classes.txt > LncRNA_upstream.txt
	
	
Downstream-lncRNA

	awk -F '\t' '{if($1==1 && $6 == "sense" && $7 == "intergenic" && $8 <= 2000 && $10 == "downstream") {print $0}}' lncRNA_classes.txt > LncRNA_downstream.txt


Intergenic-lncRNA

	awk -F '\t' '{if($1==1 && $7 == "intergenic" && $8 > 2000) {print $0}}' lncRNA_classes.txt > LncRNA_intergenic.txt


Bidirectional-lncRNA

	awk -F '\t' '{if($1==1 && $6 == "antisense" && $7 == "intergenic" && $8 <= 2000 && $10 == "upstream") {print $0}}' lncRNA_classes.txt > LncRNA_Bidirectional.txt
	
	
	
	
## **10. TE-derived lncRNAs**
cat TE.bed

	Chr1    15827287        15838845        LTR_Gypsy
	Chr1    13085455        13085593        LTR_Copia
	Chr1    11181821        11181959        LTR_Copia
	Chr1    20699111        20699248        LTR_Copia
	...
	
	
cat LncRNA.bed

	Chr1    11171031        11171031        MSTRG.2781.1
	Chr1    12199350        12199350        MSTRG.2973.1
	Chr1    13466928        13466928        MSTRG.3115.1
	Chr1    13838536        13838536        MSTRG.3127.1
	...
	
	
TE-lncRNA

	bedtools intersect -a lncRNA.bed -b TE.bed -wo | sort -u | awk '{print $1,$2,$3,$4,$6,$7,$8,$9}' | sed 's/ /\t/g' | sed '1iChr\tLncRNA_start\tLncRNA_end\tLncRNA_ID\tTE_start\tTE_end\tTE_ID\tOverlap' > TE_lncRNA_intersect.txt 
	
	
	
	
	
	
	
