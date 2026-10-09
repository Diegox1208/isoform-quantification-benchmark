#### created by Diego G. Salas
#### Test configuration, same variables as config/config.sh, pointing at one megabase of chr22 in tests/data.
#### Raw files are committed in tests/data/raw, everything derived goes to tests/results.

##### Set variables
repoDir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
dataDir="${repoDir}/tests/results"
refDir="${dataDir}/reference"
rawDir="${repoDir}/tests/data/raw"
simDir="${dataDir}/sim"
alignDir="${dataDir}/align"
quantDir="${dataDir}/quant"
resultsDir="${dataDir}"
logDir="${repoDir}/logs"
infoFile="${repoDir}/tests/samples.tsv"

##### Reference
genome="GRCh38"
chrom="chr22"
regionStart=20000001
regionEnd=21000000
fastaOffset=20000000
genomeUrl="built by tests/make_test_data.sh"
gtfUrl="built by tests/make_test_data.sh"
geneTypes="protein_coding,lncRNA"
genomeSAindexNbases=8

##### Simulation, the ground truth
truthSeed=1
fragmentsPerSample=30000
expressedGeneFraction=0.80
isoformAlpha=0.5
minTranscriptLength=400
fragmentMean=250
fragmentSd=25
errorRate=0.005

##### Benchmark
minTruth=20

##### Global parameters
threads=2
tmpDir="${TMPDIR:-/tmp}"

##### Timestamp and log
day="$(date +%Y%m%d)"
time_="$(date +%H%M%S)"
date_="${day}_${time_}"
logFile="${logDir}/test.tsv"
mkdir -p "${refDir}" "${simDir}" "${alignDir}" "${quantDir}" "${resultsDir}" "${logDir}"
