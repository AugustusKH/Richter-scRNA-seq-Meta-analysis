# This script is used to download bamfiles from ncbi using sra-tools for phs002458.v2.p1
# This dataset need request for data access before we get .ngc file for sra download

# Load sra module
module load sra-tools/3.0.3

# Data prefetch
prefetch.3 --ngc prj_40921.ngc SRR22839006
prefetch.3 --ngc prj_40921.ngc SRR22839010
prefetch.3 --ngc prj_40921.ngc SRR22839013
prefetch.3 --ngc prj_40921.ngc SRR22839018
prefetch.3 --ngc prj_40921.ngc SRR22956586
prefetch.3 --ngc prj_40921.ngc SRR22960566
prefetch.3 --ngc prj_40921.ngc SRR22960567
prefetch.3 --ngc prj_40921.ngc SRR22964505
prefetch.3 --ngc prj_40921.ngc SRR22964506
prefetch.3 --ngc prj_40921.ngc SRR22964520

# Data split files
fasterq-dump SRR22839006
fasterq-dump SRR22839010
fasterq-dump SRR22839013
fasterq-dump SRR22839018
fasterq-dump SRR22956586
fasterq-dump SRR22960566
fasterq-dump SRR22960567
fasterq-dump SRR22964505
fasterq-dump SRR22964506
fasterq-dump SRR22964520
