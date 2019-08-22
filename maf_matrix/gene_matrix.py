from __future__ import print_function
import os
from datetime import datetime
import json

import pandas as pd
import numpy as np
import requests

from .maf_object import MafObject
from .plot import MatrixPlotter, COLOR_DICT


class GeneMatrix:
    """Holds patient-gene matrix info for a single pathway.
    Write matrix to file with call to export_matrix."""

    def __init__(self, maf_obj, genes=None, annot=None, hypermutated=None,
                 n_all_patients=None, all_patients=None, trim_genes=True,
                 lookup_aa=False):
        """Builds genepatients_dict (dict): {gene: patient_list}

        Annot can be string (use column from mafObj) or dict.

        Args:
            maf_obj (MafObject)
            genes (iterable): genes to be present in matrix
            annot (dict or str): {(hugo, patient): annot (str)} or m_df column
            all_patients (list): list of all patient IDs. if not provided,
                this is deduced from maf object.
        """
        df_in = maf_obj.df
        if not all_patients:
            all_patients = sorted(list(df_in[maf_obj.patient_col].unique()))
        if not n_all_patients:
            n_all_patients = len(all_patients)

        # FILTER GENES
        df = df_in[df_in[maf_obj.hugo_col].isin(genes)].copy()
        if lookup_aa:
            df['aa_change'] = df.apply(lambda x: self.get_aa(x, maf_obj), axis=1)

        # import pdb; pdb.set_trace()
        has_class = True if maf_obj.class_col else False
        has_annot = True if annot else False
        has_aa = True if lookup_aa else False

        # BUILD GENE-SET DATAFRAME (df), inc annotation if provided
        rename_dict = {maf_obj.hugo_col: 'hugo',
                       maf_obj.patient_col: 'patient',
                       }
        use_columns = [maf_obj.hugo_col, maf_obj.patient_col]
        if has_class:
            use_columns.append(maf_obj.class_col)
            rename_dict[maf_obj.class_col] = 'var_class'
        if has_annot:
            if type(annot) == str:
                use_columns.append(annot)
                rename_dict[annot] = 'annot'
            elif type(annot) == dict:
                # left join df to dict, adding column 'annot'
                annot_series = pd.Series(annot, name='annot')
                df = df.join(annot_series, on=['patient_id', 'hugo_symbol'],
                             how='left', lsuffix='_orig')
                use_columns.append('annot')
            else:  # assume iterable of tuples ((patient, gene), annot)
                i, j = zip(*annot)
                annot_series = pd.Series(index=pd.MultiIndex.from_tuples(i),
                                         data=j)
                annot_series.name = 'annot'
                df = df.join(annot_series, on=['patient_id', 'hugo_symbol'],
                             how='left', lsuffix='_orig')
                use_columns.append('annot')
        if has_aa:
            use_columns.append('aa_change')

        df = df[use_columns].rename(columns=rename_dict)
        # ADD MISSING COLUMNS
        if not has_class:
            df['var_class'] = pd.np.nan
        if not has_annot:
            df['annot'] = pd.np.nan
        if not has_aa:
            df['aa_change'] = pd.np.nan

        # REDUCE TO ONE ROW PER PATIENT-GENE PAIR
        if has_annot or has_class:
            df = df.groupby(['patient', 'hugo']).agg({
                    'annot': self._series_to_str,
                    'var_class': lambda s: s.iloc[0],
                    'aa_change': lambda s: '|'.join([str(i) for i in list(s) if
                                                     i]),
                }).reset_index()

        df.set_index(['patient', 'hugo'], inplace=True)

        if trim_genes:
            self.genes = [i for i in genes if i in df.index.levels[1].unique()]
        else:
            self.genes = genes

        # ATTRIBUTES
        self.df = df
        self.has_class = has_class
        self.has_annot = has_annot
        # self.genes = genes
        self.genepatients = df.reset_index().groupby('hugo').apply(
            lambda s: list(s['patient'].unique())).to_dict()
        self.hypermutated = hypermutated if hypermutated is not None else []
        self.n_all_patients = n_all_patients
        self.all_patients = all_patients
        self.matrix = self._build_matrix()

    def _build_matrix(self):
        """ build gene-patient dataframe, sorted for display.
                  a      e      b      c      d      f
        KRAS      k   True      m   True   True  False
        TP53   True   True   True  False  False   True
        EGFR   True      e  False  False  False  False
        BRAF  False  False  False  False  False
        """

        patients = set()
        if self.all_patients:
            patients = self.all_patients
        else:
            for gene in self.genepatients:
                patients = patients.union(set(self.genepatients[gene]))
            patients = sorted(list(patients))

        patients = sorted(list(patients))

        # sort genes
        genes = sorted(self.genes)
        # genes.sort(key=lambda g: g in self.exclusive_genes, reverse=True)
        genes.sort(key=lambda g: len(self.genepatients[g])
            if g in self.genepatients else 0, reverse=True)

        # This is original 'matrix' input to matlab
        m = pd.DataFrame(index=genes, columns=patients, data=False)
        for gene, patient_list in self.genepatients.items():
            m.loc[gene, patient_list] = True

        # annot_dict_use = dict(
        #     [i for i in self.annot_dict.items() if i[0][0] in self.genes])
        annot_series = self.df.annot[~self.df.annot.isnull()]  # type: pd.Series
        for (patient, gene), annot in annot_series.iteritems():
            if annot:
                m.loc[gene, patient] = annot

        # m_status = ~(m == False)
        m_status = m.applymap(lambda x: x is not False)

        # RE-SORT BY WEIGHT (powers of 2)
        self.n_genes = len(self.genes)
        twos = np.power(2, range(self.n_genes - 1, -1, -1),
                        dtype=np.float64)  # gene weights

        weights = m_status.apply(lambda s: sum(s * twos), axis=0)
        patient_order = list(weights.reset_index().sort_values(
            [0, 'index'], axis=0, ascending=[False, True])['index'].values)

        matrix = m.reindex(columns=patient_order)
        return matrix

    # def export_matrix(self, outfile):
    #     """Writes tab-separated matrix file for patient/gene
    #     pairs in pathway."""
    #     pass

    def plot_matrix_orig(self, fontsize=8, max_label=150, box_px=20,
                         show_limits=False):
        mp = MatrixPlotter(self.matrix, n_all_patients=self.n_all_patients,
                           hypermutated=self.hypermutated)
        hfig, ax = mp.draw_matrix(fontsize=fontsize, max_label=max_label,
                                  show_limits=show_limits, box_px=box_px)
        return hfig, ax

    def plot_matrix(self, color_dict=None, show_title=True,
                    show_all_patients=False):
        if not color_dict:
            color_dict = COLOR_DICT
        mp = MatrixPlotter(self.matrix, self.df, color_dict=color_dict,
                           n_all_patients=self.n_all_patients,
                           show_title=show_title,
                           show_all_patients=show_all_patients)
        hfig, ax = mp.draw_matrix(fontsize=16, box_px=30)
        return hfig, ax

    def save_matrix(self, out_path=None, **plot_kwargs):
        if not out_path:
            d = datetime.utcnow()
            out_dir = os.getcwd()
            fig_name = 'matrix_' + d.strftime('%Y-%m-%d_%H%M%S%f') + '.pdf'
            out_path = os.path.join(out_dir, fig_name)
        hfig, ax = self.plot_matrix(**plot_kwargs)
        hfig.savefig(out_path)
        return out_path

    @staticmethod
    def _series_to_str(s):
        if not (~s.isnull()).any():
            return pd.np.nan
        vals = sorted(list(s.unique()))
        if len(vals) == 1:
            return str(vals[0])
        else:
            return str(vals[0]) + '+'

    def get_aa(self, d, maf_obj):
        # import pdb; pdb.set_trace()
        if not pd.isnull(d['aa_change']):
            return d['aa_change']
        else:
            return get_aa_change(chrom=d[maf_obj.chrom_col],
                                 start=d[maf_obj.start_col],
                                 end=d[maf_obj.end_col],
                                 ref=d[maf_obj.ref_col],
                                 alt=d[maf_obj.alt_col])

    def print_oncoprinter_text(self):
        """Return string formatted for oncoprint web app."""
        if not self.lookup_aa_bool:
            self.lookup_aa_bool = True
            self.df_sm = self._get_df_sm()
            self.get_matrix()
        print("\t".join(['sample', 'gene', 'alteration']))
        for gene in self.matrix.columns:
            for patient, has_mutation in iter(self.matrix[gene].items()):
                if has_mutation:
                    if self.lookup_aa_bool:
                        aa_changes = ','.join(self.df_sm.loc[
                                                  (self.df[
                                                       self.hugo_col] == gene)
                                                  & (self.df[
                                                         self.patient_col] == patient),
                                                  'aa_change'])
                        print("{}\t{}\t{}".format(patient, gene, aa_changes))
                    else:
                        print("{}\t{}\t{}".format(patient, gene, 'Y'))
        for patient in self.extra_patients:
            print(patient + '\t\t')

    def get_oncoprintjs_json(self):
        data = []
        for gene in self.matrix.columns:
            for patient, has_mutation in iter(self.matrix[gene].items()):
                d = dict(sample=patient, gene=gene)
                if has_mutation:
                    d['mutation'] = True
                data.append(d)
        return json.dumps(data)


def get_aa_change(chrom=None, start=None, end=None, ref=None, alt=None):
    vals = [chrom, start, end, ref, alt]
    if None in vals:
        raise Exception("Missing keyword argument.")
    base_url = 'http://www.broadinstitute.org/oncotator/mutation/'
    r = requests.get(base_url + '_'.join([str(i) for i in vals]))
    aa_change = r.json()['protein_change'][2:]
    return aa_change


if __name__ == "__main__" and __package__ is None:
    __package__ = "maf_matrix"
