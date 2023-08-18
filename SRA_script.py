
#!/usr/bin/python3.10
import os
import subprocess
import random
import pandas as pd

# Input parameters
df_cancer=pd.read_csv(r"/home/sam/df_control_IC.txt")
sra_names=df_cancer['Run'].tolist()
sra_reads=df_cancer['reads'].tolist()

sra_file = r"/home/sam/df_control_IC.txt"
genome_dir = 'human_ge'  # Replace with your genome directory

# Fastq-dump parameters
n_reads = 1000000  # Number of reads to download per SRA file
output_dir = '/home/sam/fastq'  # Replace with your desired output directory

# Bowtie2 parameters
num_threads = 8 
alignment_options = '--very-sensitive-local'  # Alignment options to use

# Combine BAM files using samtools
bam_files=[]
bam_file = '/home/sam/all_samples.bam'  # Replace with your desired output BAM file

bedgraph_file = '/home/sam/all_samples.bedGraph'  # Replace with your desired output bedGraph file

# Download SRA files and convert to Fastq, then align reads to the reference genome using Bowtie2

n_test=sra_names[0:7]
print(len(sra_names))
for  n in range(3,7):
    sra_name=sra_names[n]
    read=sra_reads[n]
    min_num = 0
    max_num = read
    num_samples = 1000000
    start=random.randint(0,max_num-num_samples)
    end=start+999999
    
    # Download the first n_reads from the SRA file and output as paired-end Fastq files
    output_file_1 = f'{output_dir}/{sra_name}_1.fastq.gz'
    output_file_2 = f'{output_dir}/{sra_name}_2.fastq.gz'
    bed= f'/mnt/c/users/jeremie/desktop/sra_files_IC/control_IC_{sra_name}.bed.gz'
    cmd = f'fastq-dump --gzip --split-files -N {start} -X {end} {sra_name} --outdir {output_dir} '
    subprocess.run(cmd, shell=True)

    # Align paired-end reads to the reference genome
    input_file_1 = f'{output_dir}/{sra_name}_1.fastq.gz'
    input_file_2 = f'{output_dir}/{sra_name}_2.fastq.gz'
    output_file = f'{output_dir}/{sra_name}.bam'
    cmd = f'bowtie2 -p {num_threads} -x {genome_dir} -1 {input_file_1} -2 {input_file_2} | samtools view -bS - > {output_file}'
    subprocess.run(cmd, shell=True)

    # Remove unnecessary Fastq files
    os.remove(input_file_1)
    os.remove(input_file_2)

    # Add the resulting BAM file to the list of BAM files
    bam_files.append(output_file)



# Create bedGraph file using genomeCoverageBed from bedtools
    cmd = f'samtools sort -@ {num_threads} -o {output_file} {output_file}; samtools index {output_file}; bedtools bamtobed -i {output_file} | gzip > {bed}'
    subprocess.run(cmd, shell=True)
    os.remove(output_file)
      
