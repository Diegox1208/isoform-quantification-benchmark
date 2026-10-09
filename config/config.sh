#### created by Diego G. Salas
#### Working configuration, the whole of chr22. Sourced by run.sh and by every workflow script.
#### All paths are relative to the repository root. Never edit a workflow script to change a path.

##### Set variables
repoDir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
dataDir="${repoDir}/data"
refDir="${dataDir}/reference"
rawDir="${refDir}/raw"
simDir="${dataDir}/sim"
alignDir="${dataDir}/align"
quantDir="${dataDir}/quant"
resultsDir="${repoDir}/results"
logDir="${repoDir}/logs"
infoFile="${repoDir}/config/samples.tsv"

##### Reference
genome="GRCh38"
chrom="chr22"
regionStart=""
regionEnd=""
fastaOffset=0
genomeUrl="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
gtfUrl="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz"
geneTypes="protein_coding,lncRNA"
genomeSAindexNbases=11

##### Simulation, the ground truth
truthSeed=1
fragmentsPerSample=1000000
expressedGeneFraction=0.80
isoformAlpha=0.5
minTranscriptLength=400
fragmentMean=250
fragmentSd=25
errorRate=0.005

##### Benchmark
minTruth=20

##### Global parameters
threads=6
tmpDir="${TMPDIR:-/tmp}"

##### Timestamp and log
day="$(date +%Y%m%d)"
time_="$(date +%H%M%S)"
date_="${day}_${time_}"
logFile="${logDir}/run.tsv"
mkdir -p "${refDir}" "${simDir}" "${alignDir}" "${quantDir}" "${resultsDir}" "${logDir}"
