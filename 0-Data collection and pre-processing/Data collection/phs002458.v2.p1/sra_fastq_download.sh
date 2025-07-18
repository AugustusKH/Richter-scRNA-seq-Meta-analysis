# This script is used to download bamfiles from ncbi using sra-tools for phs002458.v2.p1
# This dataset need request for data access before we get .ngc file for sra download

# Load sra module
module load sra-tools/3.0.3

# Data prefetch
prefetch.3 --ngc prj_40921.ngc SRR22882959
prefetch.3 --ngc prj_40921.ngc SRR22882960
prefetch.3 --ngc prj_40921.ngc SRR22882961
prefetch.3 --ngc prj_40921.ngc SRR22882962
prefetch.3 --ngc prj_40921.ngc SRR22882969
prefetch.3 --ngc prj_40921.ngc SRR22882970
prefetch.3 --ngc prj_40921.ngc SRR22882971
prefetch.3 --ngc prj_40921.ngc SRR22882972
prefetch.3 --ngc prj_40921.ngc SRR22882973
prefetch.3 --ngc prj_40921.ngc SRR22882974
prefetch.3 --ngc prj_40921.ngc SRR22882975
prefetch.3 --ngc prj_40921.ngc SRR22882976
prefetch.3 --ngc prj_40921.ngc SRR22882977
prefetch.3 --ngc prj_40921.ngc SRR22882978
prefetch.3 --ngc prj_40921.ngc SRR22882979
prefetch.3 --ngc prj_40921.ngc SRR22882980
prefetch.3 --ngc prj_40921.ngc SRR22882981
prefetch.3 --ngc prj_40921.ngc SRR22882982
prefetch.3 --ngc prj_40921.ngc SRR22882984
prefetch.3 --ngc prj_40921.ngc SRR22883041
prefetch.3 --ngc prj_40921.ngc SRR22883042
prefetch.3 --ngc prj_40921.ngc SRR22883043
prefetch.3 --ngc prj_40921.ngc SRR22883044
prefetch.3 --ngc prj_40921.ngc SRR22883045
prefetch.3 --ngc prj_40921.ngc SRR22883046
prefetch.3 --ngc prj_40921.ngc SRR22883047
prefetch.3 --ngc prj_40921.ngc SRR22883048
prefetch.3 --ngc prj_40921.ngc SRR22883052
prefetch.3 --ngc prj_40921.ngc SRR22883053
prefetch.3 --ngc prj_40921.ngc SRR22883054
prefetch.3 --ngc prj_40921.ngc SRR22883055
prefetch.3 --ngc prj_40921.ngc SRR22883056
prefetch.3 --ngc prj_40921.ngc SRR22883057
prefetch.3 --ngc prj_40921.ngc SRR22883058
prefetch.3 --ngc prj_40921.ngc SRR22883059
prefetch.3 --ngc prj_40921.ngc SRR22883060
prefetch.3 --ngc prj_40921.ngc SRR22883061
prefetch.3 --ngc prj_40921.ngc SRR22903814
prefetch.3 --ngc prj_40921.ngc SRR22903815

# Data split files
fasterq-dump SRR22882959
fasterq-dump SRR22882960
fasterq-dump SRR22882961
fasterq-dump SRR22882962
fasterq-dump SRR22882969
fasterq-dump SRR22882970
fasterq-dump SRR22882971
fasterq-dump SRR22882972
fasterq-dump SRR22882973
fasterq-dump SRR22882974
fasterq-dump SRR22882975
fasterq-dump SRR22882976
fasterq-dump SRR22882977
fasterq-dump SRR22882978
fasterq-dump SRR22882979
fasterq-dump SRR22882980
fasterq-dump SRR22882981
fasterq-dump SRR22882982
fasterq-dump SRR22882984
fasterq-dump SRR22883041
fasterq-dump SRR22883042
fasterq-dump SRR22883043
fasterq-dump SRR22883044
fasterq-dump SRR22883045
fasterq-dump SRR22883046
fasterq-dump SRR22883047
fasterq-dump SRR22883048
fasterq-dump SRR22883052
fasterq-dump SRR22883053
fasterq-dump SRR22883054
fasterq-dump SRR22883055
fasterq-dump SRR22883056
fasterq-dump SRR22883057
fasterq-dump SRR22883058
fasterq-dump SRR22883059
fasterq-dump SRR22883060
fasterq-dump SRR22883061
fasterq-dump SRR22903814
fasterq-dump SRR22903815







