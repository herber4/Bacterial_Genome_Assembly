mamba create -c bioconda -c conda-forge -n trycycler trycycler
mamba install bioconda::flye
mamba install bioconda::raven-assembler
pip install medaka
pip uninstall ont-fast5-api
pip install ont-fast5-api


### trycycler subsample reads

trycycler subsample --reads reads.fastq --genome_size 5.2m --threads 8 --out_dir read_subsets


# here is draft genome assembly using Raven and FLYE

#!/bin/bash
#SBATCH --job-name=assemble
#SBATCH -n 24
#SBATCH --partition=bigmem
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --output=/path/to/assembly/data/assemb.%j.out
#SBATCH --error=/path/to/assembly/data/assemb.%j.err
#SBATCH --mail-user=user@clemson.edu
# Load change dir
cd /path/to/assembly/data/

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh

conda activate trycycler

threads=24  # change as appropriate for your system
mkdir assemblies

export PATH=/data2/lackey_lab/austin/peach/new_copper/tri:$PATH

flye --nano-hq read_subsets/sample_01.fastq --threads "$threads" --out-dir assembly_01 && cp assembly_01/assembly.fasta assemblies/assembly_01.fasta && cp assembly_01/assembly_graph.gfa assemblies/assembly_01.gfa && rm -r assembly_01
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_03.gfa read_subsets/sample_03.fastq > assemblies/assembly_03.fasta

flye --nano-hq read_subsets/sample_04.fastq --threads "$threads" --out-dir assembly_04 && cp assembly_04/assembly.fasta assemblies/assembly_04.fasta && cp assembly_04/assembly_graph.gfa assemblies/assembly_04.gfa && rm -r assembly_04
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_06.gfa read_subsets/sample_06.fastq > assemblies/assembly_06.fasta

flye --nano-hq read_subsets/sample_07.fastq --threads "$threads" --out-dir assembly_07 && cp assembly_07/assembly.fasta assemblies/assembly_07.fasta && cp assembly_07/assembly_graph.gfa assemblies/assembly_07.gfa && rm -r assembly_07
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_09.gfa read_subsets/sample_09.fastq > assemblies/assembly_09.fasta

flye --nano-hq read_subsets/sample_10.fastq --threads "$threads" --out-dir assembly_10 && cp assembly_10/assembly.fasta assemblies/assembly_10.fasta && cp assembly_10/assembly_graph.gfa assemblies/assembly_10.gfa && rm -r assembly_10
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_12.gfa read_subsets/sample_12.fastq > assemblies/assembly_12.fasta


# here is clustering step of trycycler

#!/bin/bash
#SBATCH --job-name=cluster
#SBATCH -n 24
#SBATCH --partition=bigmem
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --output=/path/to/assembly/data/cluster.%j.out
#SBATCH --error=/path/to/assembly/data/cluster.%j.err
#SBATCH --mail-user=user@clemson.edu
# Load change dir
cd /path/to/assembly/data/

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh

conda activate trycycler

trycycler cluster --assemblies assemblies/*.fasta --reads reads.fastq --out_dir trycycler --threads 24


# here is reconcile step of trycycler

#!/bin/bash
#SBATCH --job-name=reconcile
#SBATCH -n 24
#SBATCH --partition=bigmem
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --output=/path/to/assembly/data/reconcile.%j.out
#SBATCH --error=/path/to/assembly/data/reconcile.%j.err
#SBATCH --mail-user=user@clemson.edu
# Load change dir
cd /path/to/assembly/data/

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh

conda activate trycycler

trycycler reconcile --threads 24 --min_1kbp_identity 23 --reads reads.fastq --cluster_dir trycycler/cluster_001
trycycler reconcile --threads 24 --min_1kbp_identity 23 --reads reads.fastq --cluster_dir trycycler/cluster_002
trycycler reconcile --threads 24 --min_1kbp_identity 23 --reads reads.fastq --cluster_dir trycycler/cluster_003


# here is MSA step of trycycler

#!/bin/bash
#SBATCH --job-name=msa
#SBATCH -n 24
#SBATCH --partition=bigmem
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --output=/path/to/assembly/data/msa.%j.out
#SBATCH --error=/path/to/assembly/data/msa.%j.err
#SBATCH --mail-user=user@clemson.edu
# Load change dir
cd /path/to/assembly/data/

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh

conda activate trycycler

trycycler msa --threads 24 --cluster_dir trycycler/cluster_001
trycycler msa --threads 24 --cluster_dir trycycler/cluster_002
trycycler msa --threads 24 --cluster_dir trycycler/cluster_003


# here is the partition step of trycycler

#!/bin/bash
#SBATCH --job-name=partition
#SBATCH -n 24
#SBATCH --partition=bigmem
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --output=/path/to/assembly/data/partition.%j.out
#SBATCH --error=/path/to/assembly/data/partition.%j.err
#SBATCH --mail-user=user@clemson.edu
# Load change dir
cd /path/to/assembly/data/

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh

conda activate trycycler

trycycler partition --threads 24 --reads reads.fastq --cluster_dirs trycycler/cluster_001 trycycler/cluster_002 trycycler/cluster_003


# here is the consensus step of strycycler

#!/bin/bash
#SBATCH --job-name=consensus
#SBATCH -n 24
#SBATCH --partition=bigmem
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --output=/path/to/assembly/data/consensus.%j.out
#SBATCH --error=/path/to/assembly/data/consensus.%j.err
#SBATCH --mail-user=user@clemson.edu
# Load change dir
cd /path/to/assembly/data/

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh

conda activate trycycler
trycycler consensus --threads 24 --cluster_dir trycycler/cluster_001
trycycler consensus --threads 24 --cluster_dir trycycler/cluster_002
trycycler consensus --threads 24 --cluster_dir trycycler/cluster_003


cat trycycler/cluster_*/7_final_consensus.fasta > trycycler/consensus.fasta


medaka_consensus -i 4_reads.fastq -d 7_final_consensus.fasta -o ./medaka -t 8 --bacteria
medaka_consensus -i 4_reads.fastq -d 7_final_consensus.fasta -o ./medaka -t 8 --bacteria
medaka_consensus -i 4_reads.fastq -d 7_final_consensus.fasta -o ./medaka -t 8 --bacteria

cat cluster_*/medaka/medaka_consensus > medaka_consensus.fasta


# here i will do polypolish
mamba install bioconda::polypolish

ml bwa

bwa index medaka_consensus.fasta
bwa mem -t 8 -a medaka_consensus.fasta CuR_R1_clean.fastq.gz > alignments_1.sam
bwa mem -t 8 -a medaka_consensus.fasta CuR_R2_clean.fastq.gz > alignments_2.sam
polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish polish medaka_consensus.fasta filtered_1.sam filtered_2.sam > polypolished_consensus.fasta


# now to do meryl and merqury QC assessment
meryl k=16 count output r1.meryl R1_R1.fastq.gz
meryl k=16 count output r2.meryl R1_R2.fastq.gz

#merge meryl db

meryl union-sum output union.meryl r1.meryl r2.meryl/

#now run merqury

../merqury/./merqury.sh union.meryl R1_COMPLETE_GENOME.fasta testr
merqury.sh union.meryl polypolished_consensus.fasta test


# busco analysis
conda create --name busco python=3.7
conda activate busco
conda install -c bioconda conda-forge -c bioconda busco=5.2.2

conda activate busco

busco --download prokaryota

busco --offline -c 8 -l busco_downloads/lineages/xanthomonadales_odb10/ -i polypolished_consensus.fasta -o busco_out


# prokka annotation of genes and proteins

prokka --outdir results --prefix Xap_CuR --locustag --centre --compliant polypolished_consensus.fasta

