"""
Copyright (C) 2015 Stephen Gaffney

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
__author__ = 'Stephen G. Gaffney'

import os
import json
import argparse
from collections import defaultdict
from datetime import datetime

import pandas as pd
import requests
import numpy as np

from maf_matrix.plot import MatrixPlotter

test_maf_path = 'chol_combined_removed_hypermutated_AA39.maf'
test_use_genes = ['KRAS', 'TP53', 'MUC4']
test_lookup_aa_bool = False


def get_aa_change(chrom=None, start=None, end=None, ref=None, alt=None):
    vals = [chrom, start, end, ref, alt]
    if None in vals:
        raise Exception("Missing keyword argument.")
    base_url = 'http://www.broadinstitute.org/oncotator/mutation/'
    r = requests.get(base_url + '_'.join([str(i) for i in vals]))
    aa_change = r.json()['protein_change'][2:]
    return aa_change


def get_patient_from_barcode(barcode):
    if barcode.startswith('TCGA'):
        return barcode[:12]
    else:
        return barcode


class MafObject(object):
    """Holds maf as dataframe.

    Dataframe df is set once.
    Dataframe df_sm is re-set with each update to use_genes
    """
    def __init__(self, maf_path_or_df, lookup_aa_bool=False, show_all_patients=False):
        # TEST ARGUMENTS
        if type(maf_path_or_df) == str:
            if not os.path.isfile(maf_path_or_df):
                raise ValueError("Invalid maf path.")
            self.maf_path = maf_path_or_df
            df_in = None
        elif type(maf_path_or_df) == pd.core.frame.DataFrame:
            df_in = maf_path_or_df

        self.lookup_aa_bool = lookup_aa_bool
        self.show_all_patients = show_all_patients
        self.patient_col = None
        self.hugo_col = None
        self.chrom_col = None
        self.class_col = None
        self.ref_col = None
        self.alt_col = None
        self.start_col = None
        self.end_col = None
        self.n_all_patients = None
        self.df = self._get_maf_df(df_in)

        # gene-list dependent
        self._use_genes = None
        self.use_patients = None
        self.matrix = None
        self.df_sm = None
        self.extra_patients = None

    @property
    def use_genes(self):
        return self._use_genes

    @use_genes.setter
    def use_genes(self, genes):
        """Takes iterator object with hugo symbols."""
        self._use_genes = genes
        self.df_sm = self._get_df_sm()
        self.matrix = self.get_matrix()

    def _get_maf_df(self, df_in=None):
        """Populate column name attributes and whole-maf dataframe."""
        if df_in is None:
            read_kwargs = dict(sep='\t', comment='#')
            df = pd.read_csv(self.maf_path, **read_kwargs)
        else:
            df = df_in
        df.rename(columns=lambda c: c.lower(), inplace=True)
        # GET COLUMN NAMES
        hugo_col = [i for i in df.columns if 'hugo' in i][0]
        chrom_col = [i for i in df.columns if 'chrom' in i][0]
        class_col = [i for i in df.columns if 'classification' in i][0]
        ref_col = [i for i in df.columns if 'reference' in i][0]
        start_col = [i for i in df.columns if 'start' in i][0]
        end_col = [i for i in df.columns if 'end' in i][0]
        alt_col = 'alt'

        alt_cols = [i for i in df.columns if 'tumor_seq' in i]
        df[alt_col] = df.apply(lambda d: [n for n in d[alt_cols]
                                          if n != d[ref_col]][0],
                               axis=1)
        has_patient = [i for i in df.columns if 'patient' in i]
        if has_patient:
            patient_col = has_patient[0]
        else:
            barcode_col = 'tumor_sample_barcode'
            patient_col = 'patient_id'
            df[patient_col] = df[barcode_col].apply(get_patient_from_barcode).\
                astype(str)
        df[class_col] = df[class_col].str.lower()
        df[[hugo_col, patient_col, class_col]].head()
        df['aa_change'] = np.nan
        self.patient_col = patient_col
        self.hugo_col = hugo_col
        self.chrom_col = chrom_col
        self.class_col = class_col
        self.ref_col = ref_col
        self.alt_col = alt_col
        self.start_col = start_col
        self.end_col = end_col
        return df

    def _get_df_sm(self):
        use_cols = [self.hugo_col, self.patient_col, self.class_col,
                    self.chrom_col, self.start_col, self.end_col,
                    self.ref_col, self.alt_col, 'aa_change']
        df_sm = self.df.loc[
            (~self.df[self.class_col].isin(['silent', 'coding-silent']))
            & (self.df[self.hugo_col].isin(self.use_genes)), use_cols]

        if self.lookup_aa_bool:
            # ADD AA_CHANGE COLUMN
            df_sm['aa_change'] = df_sm.apply(
                lambda d: self.__get_aa_change_for_row(d), axis=1)
            # COPY TO DF
            for ind, row in df_sm.iterrows():
                self.df.loc[ind, 'aa_change'] = row.aa_change
        return df_sm

    def __get_aa_change_for_row(self, d):
        if not pd.isnull(d['aa_change']):
            return d['aa_change']
        else:
            return get_aa_change(chrom=d[self.chrom_col],
                                 start=d[self.start_col],
                                 end=d[self.end_col],
                                 ref=d[self.ref_col],
                                 alt=d[self.alt_col])

    def get_matrix(self):

        self.use_genes.sort()
        use_patients = list(self.df_sm['patient_id'].unique())
        use_patients.sort()
        n_patients_dict = defaultdict(
            int, self.df_sm.groupby(self.hugo_col)[self.patient_col].nunique())
        n_genes_dict = dict(
            self.df_sm.groupby(self.patient_col)[self.hugo_col].nunique())
        self.use_genes.sort(key=lambda j: n_patients_dict[j], reverse=True)
        use_patients.sort(key=lambda k: n_genes_dict[k], reverse=True)
        matrix = pd.DataFrame(index=use_patients, columns=self.use_genes,
                              dtype=bool, data=False)
        for ind, row in self.df_sm.iterrows():
            matrix.loc[row.patient_id, row.hugo_symbol] = True

        # RE-SORT BY WEIGHT (powers of 2)
        n_genes = len(self.use_genes)
        twos = np.power(2, xrange(n_genes-1, -1, -1))  # gene weights
        twos_dict = dict(matrix.apply(lambda r: sum(r*twos), axis=1))
        use_patients.sort(key=lambda i: twos_dict[i], reverse=True)
        self.use_patients = use_patients
        self.extra_patients = list(set(self.df[self.patient_col].unique())
                                   .difference(set(self.use_patients)))
        self.extra_patients.sort()
        # REINDEX USING WEIGHTS
        if self.show_all_patients:
            matrix = matrix.reindex(self.use_patients + self.extra_patients,
                                    fill_value=False)
        else:
            matrix = matrix.reindex(self.use_patients)
        self.n_all_patients = len(self.use_patients) + len(self.extra_patients)
        return matrix

    def print_oncoprinter_text(self):
        """Return string formatted for oncoprint web app."""
        if not self.lookup_aa_bool:
            self.lookup_aa_bool = True
            self.df_sm = self._get_df_sm()
            self.get_matrix()
        print "\t".join(['sample', 'gene', 'alteration'])
        for gene in self.matrix.columns:
            for patient, has_mutation in self.matrix[gene].iteritems():
                if has_mutation:
                    if self.lookup_aa_bool:
                        aa_changes = ','.join(self.df_sm.loc[
                            (self.df[self.hugo_col] == gene)
                            & (self.df[self.patient_col] == patient),
                            'aa_change'])
                        print "{}\t{}\t{}".format(patient, gene, aa_changes)
                    else:
                        print "{}\t{}\t{}".format(patient, gene, 'Y')
        for patient in self.extra_patients:
            print patient + '\t\t'

    def get_oncoprintjs_json(self):
        data = []
        for gene in self.matrix.columns:
            for patient, has_mutation in self.matrix[gene].iteritems():
                d = dict(sample=patient, gene=gene)
                if has_mutation:
                    d['mutation'] = True
                data.append(d)
        return json.dumps(data)

    def plot_matrix(self):
        mp = MatrixPlotter(self.matrix, n_all_patients=self.n_all_patients)
        fig1 = mp.draw_matrix()
        return fig1

    def save_matrix(self, out_path=None):
        if not out_path:
            d = datetime.utcnow()
            maf_dir = os.path.dirname(self.maf_path) \
                if self.maf_path else os.getcwd()
            fig_name = 'matrix_' + d.strftime('%Y-%m-%d_%H%M%S%f') + '.pdf'
            out_path = os.path.join(maf_dir, fig_name)
        fig1 = self.plot_matrix()
        fig1.savefig(out_path, bbox_inches='tight')
        return out_path


if __name__ == '__main__':

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
        m = MafObject(test_maf_path, lookup_aa_bool=True)
        m.use_genes = args.genes
        print("Paste at http://www.cbioportal.org/oncoprinter.jsp")
        m.print_oncoprinter_text()
    else:
        if outfile_exists:
            overwrite = raw_input("Overwrite {}? (y/n): ".format(args.out_path))
            if overwrite.lower() == 'y':
                overwrite_bool = True

        if not outfile_exists or overwrite_bool:
            m = MafObject(test_maf_path, lookup_aa_bool=False,
                          show_all_patients=args.all_patients)
            m.use_genes = args.genes
            fig_path = m.save_matrix(args.out_path)
            print("Saved matrix to {}".format(fig_path))
        else:
            print("Matrix save aborted.")
