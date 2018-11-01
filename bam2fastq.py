"""This script aims to extract Illumina read pair reads from a BAM file into
separate fastq files. Specifically it aims to handle cases where reads have
been duplicated for whatever reason.
"""
import sys
import argparse
import logging
from pathlib import Path
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    import pysam


def main():
    args = cli()

    logging.info("Creating and loading index into memory...")
    index = index_bam(args.input)
    logging.info("Index loaded.")

    logging.info("Geting set of unique read ids...")
    with pysam.AlignmentFile(args.input) as bam:
        read_ids = get_read_ids(bam)
    logging.info("Gathered all unique read ids.")

    # extract unique read pairs and write to file
    read1_filepath = args.output / (args.prefix + '_r1.fastq')
    read1_fh = open(read1_filepath, 'w')
    read2_filepath = args.output / (args.prefix + '_r2.fastq')
    read2_fh = open(read2_filepath, 'w')

    logging.info("Writing unique reads to file...")
    for i, read_id in enumerate(read_ids):
        reads = index.find(read_id)
        r1s, r2s = get_unique_reads_pairs(reads)

        write_reads_to_fastq(r1s, read1_fh)
        write_reads_to_fastq(r2s, read2_fh)
        if not args.no_progress_bar:
            update_progress(i / len(read_ids))

    if not args.no_progress_bar:
        update_progress(1)

    read1_fh.close()
    read2_fh.close()


def write_reads_to_fastq(records, fh):
    if not records:
        return

    strings = []
    for r in records:
        fastq_string = to_fastq_string(r)
        if fastq_string:
            strings.append(fastq_string)

    fastq_records = "\n".join(strings)
    print(fastq_records, file=fh)


def to_fastq_string(record):
    """Extracts the information required for a fastq entry from the BAM entry

    :param record: BAM record
    :return: A string representation of the four fastq lines
    """
    qual = record.query_qualities
    seq = record.query_sequence

    if seq is None or qual is None:
        return ''
    else:
        qual = ''.join([chr(q + 33) for q in record.query_qualities])

    rg = ''
    if record.is_read1:
        rg = '/1'
    elif record.is_read2:
        rg = '/2'

    try:
        tag = f"RG:{record.get_tag('RG')}"
    except KeyError:
        tag = ''

    header = f"@{record.query_name}{rg} {tag}"

    assert len(seq) == len(qual)
    return f"{header}\n{seq}\n+\n{qual}"


def get_unique_reads_pairs(reads):
    """Finds all unique read group paired reads.

    :param reads: An iterable object containig pysam AlignmentRead objects.
    :return: Tuple of (unique read1s, unique read2s)
    """
    seen_read1_groups = set()
    seen_read2_groups = set()
    r1s = []
    r2s = []
    for read in reads:
        tag = read.get_tag('RG')
        is_unique_read1 = (read.is_read1 and tag not in seen_read1_groups)
        is_unique_read2 = (read.is_read2 and tag not in seen_read2_groups)
        if is_unique_read1:
            r1s.append(read)
            seen_read1_groups.add(tag)
        elif is_unique_read2:
            r2s.append(read)
            seen_read2_groups.add(tag)

    assert len(r1s) == len(r2s), \
        f"Missing one member of read pair for {read.query_name}"

    return r1s, r2s


def index_bam(filepath):
    """Creates an in-memory index for BAM file.

    :param filepath: Path to the BAM file.
    """
    index = pysam.IndexedReads(pysam.AlignmentFile(filepath))
    index.build()
    return index


def get_read_ids(bam):
    """Gets all unique read ids in a bam file.

    :param bam: pysam file handle for bam file
    :return: Set of read ids as strings
    """
    read_ids = set()
    for record in bam:
        if record.is_paired:
            read_ids.add(record.query_name)

    return read_ids


def handle_file(parser, arg):
    """Ensures file passed as input is a file that exists.

    :param parser: argparse parser object
    :param arg: argument passed as input
    :return: The filepath as a Path object
    """
    p = Path(arg)
    if not p.is_file():
        parser.error(f"Input file '{arg}' is not a file.")
    else:
        return p


def handle_output_dir(parser, arg):
    """Ensures the output directory exists.

    :param parser:  argparse parser object
    :param arg: argument passed as output
    :return: The output directory path as a Path object
    """
    p = Path(arg)
    if p.is_dir():
        return p
    else:
        parser.error(f"Output directory '{arg}' does not exist.")


def setup_logging(level):
    """Sets up the logging based on cli option passed."""
    logging_levels = {
        0: "NOTSET",
        1: "CRITICAL",
        2: "ERROR",
        3: "WARNING",
        4: "INFO",
        5: "DEBUG"
    }
    log_level = logging_levels.get(level)
    logging.basicConfig(level=log_level,
                        format='[%(asctime)s]:%(levelname)s:%(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')


def update_progress(progress, write_progress_bar_to=sys.stdout):
    """Creates and updates a progress bar.
    Recognition to https://stackoverflow.com/a/15860757/5299417
    :param progress: Value between 0 and 1 (percent as decimal)
    :param write_progress_bar_to: Where to write progress bar to.
    """
    bar_length = 40  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(bar_length * progress))
    progress_percent = round(progress * 100, 2)
    text = "\rPercent: [{0}] {1}% {2}".format(
        "#" * block + "-" * (bar_length - block), progress_percent, status)
    write_progress_bar_to.write(text)
    write_progress_bar_to.flush()


def cli():
    """Create command line interface and parse arguments for main program."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "-i", "--input",
        type=lambda arg: handle_file(parser, arg),
        help="""Path to BAM file.""",
        required=True)

    parser.add_argument(
        "-o", "--output",
        type=lambda arg: handle_output_dir(parser, arg),
        help="""Directory to write read pair fastqs to. Default is the 
        current directory.""",
        default=Path().cwd())

    parser.add_argument(
        "-p", "--prefix",
        type=str,
        help="""Path and prefix to write read pair fastq to. Default is the 
            basename of the BAM file. _r1 and _r2 will be added to this 
            prefix.""",
        default=None)

    parser.add_argument(
        "--log_level",
        help="Level of logging. 0 is none, 5 is for debugging. Default is 4 "
             "which will report info, warnings, errors, and critical "
             "information.",
        default=4,
        type=int,
        choices=range(6))

    parser.add_argument(
        "--no_progress_bar",
        help="Do not display progress bar.",
        action='store_true',
        default=False
    )

    args = parser.parse_args()
    setup_logging(args.log_level)

    if args.prefix is None:
        args.prefix = args.input.stem

    logging.debug(f"--input = {args.input}")
    logging.debug(f"--output = {args.output}")
    logging.debug(f"--prefix = {args.prefix}")
    logging.debug(f"--log_level = {args.log_level}")
    logging.debug(f"--no_progress_bar = {args.no_progress_bar}")

    return args


if __name__ == '__main__':
    main()
