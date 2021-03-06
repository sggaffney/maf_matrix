from .gene_matrix import MafObject, GeneMatrix


__author__ = 'Stephen G. Gaffney'


def run():
    import os
    import argparse
    # noinspection PyCompatibility
    from builtins import input

    def extant_file(file_path):
        """'Type' for argparse - checks that file exists but does not open."""
        if not os.path.isfile(file_path):
            raise argparse.ArgumentTypeError(
                "{0} does not exist".format(file_path))
        return file_path

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--maf_path",
                        help="Path to maf file.",
                        type=extant_file, required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-t', '--oncoprinter', action='store_true')
    group.add_argument("-o", "--out_path",
                       help="Output figure path including extension.",
                       type=str, default=None)
    parser.add_argument("-a", "--all_patients",
                        help="Show all patients in matrix plot.",
                        action='store_true')
    parser.add_argument("-g", "--genes",
                        help="Gene hugo symbols, case sensitive.",
                        type=str, nargs='+')
    args = parser.parse_args()

    outfile_exists = os.path.exists(args.out_path) if args.out_path else False
    overwrite_bool = False

    if args.oncoprinter:
        m = MafObject(args.maf_path, lookup_aa_bool=True)
        m.use_genes = args.genes
        print("Paste at http://www.cbioportal.org/oncoprinter.jsp")
        m.print_oncoprinter_text()
    else:
        if outfile_exists:
            overwrite = input("Overwrite {}? (y/n): ".format(args.out_path))
            if overwrite.lower() == 'y':
                overwrite_bool = True

        if not outfile_exists or overwrite_bool:
            m = MafObject(args.maf_path, lookup_aa_bool=False,
                          show_all_patients=args.all_patients)
            m.use_genes = args.genes
            fig_path = m.save_matrix(args.out_path)
            print("Saved matrix to {}".format(fig_path))
        else:
            print("Matrix save aborted.")


if __name__ == '__main__':
    pass
