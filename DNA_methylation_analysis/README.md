These scripts are intended for bisulfite sequencing data.

The first script (to be run on command line) runs demultiplexing of samples (with fumi tool), trimming & fastQC (with trim galore, which also runs cut adapt), alignment (with Bismark) and PCR deduplication (with fumi tool).

The second script (to be run on R) conducts methylation calling and DMR ID with methylkit pacakge.
