#!/usr/bin/env sh

#Rscript /orozco/services/Rdata/Web/apps/nucleServ/bin/readBAM.R --type paired --input /orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/uploads/rep3_00m_G1.bam --output /orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/uploads/rep3_00m_G1.bam.RData >> /orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/.tmp/PP_rep3_00m_G1.bam.log 2>&1

script="/orozco/services/Rdata/Web/apps/nucleServ/bin/readBAM.R"
type="paired"
input="/orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/uploads/rep3_00m_G1.bam"
output="/home/rilla/scratch/tmp/rep3_00m_G1.bam.RData"

cmd="Rscript $script --type $type --input $input --output $output"

/usr/bin/time -v $cmd

#User time (seconds): 3503.32
#System time (seconds): 17.66
#Percent of CPU this job got: 99%
#Elapsed (wall clock) time (h:mm:ss or m:ss): 59:04.13
#Average shared text size (kbytes): 0
#Average unshared data size (kbytes): 0
#Average stack size (kbytes): 0
#Average total size (kbytes): 0
#Maximum resident set size (kbytes): 51006336
#Average resident set size (kbytes): 0
#Major (requiring I/O) page faults: 0
#Minor (reclaiming a frame) page faults: 1218654
#Voluntary context switches: 2489
#Involuntary context switches: 12419
#Swaps: 0
#File system inputs: 6016376
#File system outputs: 42432
#Socket messages sent: 0
#Socket messages received: 0
#Signals delivered: 0
#Page size (bytes): 4096
#Exit status: 0

