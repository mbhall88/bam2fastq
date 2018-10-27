# bam2fastq

This script if for extracting Illumina read pairs from a BAM file into separate fastq files for each pair.  

Specifically it handles duplications of reads (even if they're not marked as duplicates) and will filter all these out before writing to file.  

If there are multiple read groups (`RG`) within the file this is also accounted for. For example, if two reads happen to have the same read ID, but are from different read groups, then they will *both* be written to file.  

## Usage

Script is compatible with Python 3.5+  

Install the single dependency `pysam`

```sh
pip3 install pysam
```

Get the script (and test script) and run tests
```sh
git clone https://github.com/mbhall88/bam2fastq.git
cd bam2fastq
python3 test_bam2fastq.py
python3 bam2fastq.py -h
```

```
usage: bam2fastq.py [-h] -i INPUT [-o OUTPUT] [-p PREFIX]
                    [--log_level {0,1,2,3,4,5}] [--no_progress_bar]

This script aims to extract Illumina read pair reads from a BAM file into
separate fastq files. Specifically it aims to handle cases where reads have
been duplicated for whatever reason.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to BAM file.
  -o OUTPUT, --output OUTPUT
                        Directory to write read pair fastqs to. Default is the
                        current directory.
  -p PREFIX, --prefix PREFIX
                        Path and prefix to write read pair fastq to. Default
                        is the basename of the BAM file. _r1 and _r2 will be
                        added to this prefix.
  --log_level {0,1,2,3,4,5}
                        Level of logging. 0 is none, 5 is for debugging.
                        Default is 4 which will report info, warnings, errors,
                        and critical information.
  --no_progress_bar     Do not display progress bar.
```
