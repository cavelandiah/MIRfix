def parseargs():
    parser = argparse.ArgumentParser(description='MIRfix automatically curates miRNA datasets by improving alignments of their precursors, the consistency of the annotation of mature miR and miR* sequence, and the phylogenetic coverage. MIRfix produces alignments that are comparable across families and sets the stage for improved homology search as well as quantitative analyses.')
    parser.add_argument("-j", "--cores", type=int, default=1, help='Number of parallel processes to run')
    parser.add_argument("-f", "--families", type=str, required=True, help='Path to list of families to work on')
    #parser.add_argument("-i", "--famdir", type=str, required=True, help='Directory where family files are located')
    parser.add_argument("-g", "--genomes", type=str, required=True, help='Genome FASTA files to parse')
    #parser.add_argument("-m", "--mapping", type=str, required=True, help='Mapping between precursor and families')
    parser.add_argument("-a", "--mature", type=str, required=True, help='FASTA files containing mature sequences')
    #parser.add_argument("-d", "--maturedir", type=str, default='', help='Directory of matures')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory for output')
    parser.add_argument("--force", action='store_true', help='Force MIRfix to overwrite existing output directories')
    parser.add_argument("-e", "--extension", type=int, default=10, help='Extension of nucleotides for precursor cutting')
    parser.add_argument("-l", "--logdir", type=str, default='LOGS', help='Directory to write logfiles to')
    parser.add_argument("--loglevel", type=str, default='WARNING', choices=['WARNING','ERROR','INFO','DEBUG','CRITICAL'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

