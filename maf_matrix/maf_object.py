import os
from collections import OrderedDict

import pandas as pd
import numpy as np


class MafObject(object):
    """Holds maf as dataframe.

    Dataframe df is set once.
    Dataframe df_sm is re-set with each update to use_genes
    """
    def __init__(self, maf_path_or_df):
        # TEST ARGUMENTS
        df_in = None
        if type(maf_path_or_df) == str:
            if not os.path.isfile(maf_path_or_df):
                raise ValueError("Invalid maf path.")
            self.maf_path = maf_path_or_df
        elif type(maf_path_or_df) == pd.DataFrame:
            df_in = maf_path_or_df

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

    def _get_maf_df(self, df_in=None):
        """Populate column name attributes and whole-maf dataframe."""
        if df_in is None:
            read_kwargs = dict(sep='\t', comment='#')
            df = pd.read_csv(self.maf_path, **read_kwargs)
        else:
            df = df_in.copy()
        df.rename(columns=lambda c: c.lower(), inplace=True)
        # GET COLUMN NAMES
        hugo_col = [i for i in df.columns if 'hugo' in i
                    or 'symbol' in i or 'gene' in i][0]
        chrom_col = [i for i in df.columns if 'chrom' in i
                     or i.startswith('chr')][0]
        class_col = [i for i in df.columns if 'classification' in i][0]
        ref_col = [i for i in df.columns if 'reference' in i or
                   i.startswith('ref')][0]
        start_col = [i for i in df.columns if 'start' in i]
        if not start_col:
            start_col = [i for i in df.columns if i.startswith('pos')][0]
        else:
            start_col = start_col[0]
        end_col = [i for i in df.columns if 'end' in i]
        end_col = end_col[0] if end_col else start_col

        alt_cols = [i for i in df.columns if 'tumor_seq' in i]
        sample_col = [i for i in df.columns if 'tumor_sample_barcode' in i or
                      i.startswith('sample')]
        sample_col = sample_col[0] if sample_col else None
        patient_col = [i for i in df.columns if i.startswith('patient')
                       or i.startswith('case')]

        if alt_cols:
            alt_col = 'alt'
            df['alt'] = df.apply(lambda d: [n for n in d[alt_cols]
                                              if n != d[ref_col]][0],
                                   axis=1)
        else:
            alt_col = [i for i in df.columns if i.startswith('alternative') or
                       i.startswith('alt')][0]

        if patient_col:
            patient_col = patient_col[0]
        else:
            barcode_col = 'tumor_sample_barcode'
            patient_col = 'patient_id'
            df[patient_col] = df[barcode_col].apply(get_patient_from_barcode).\
                astype(str)
        df[class_col] = df[class_col].str.lower()
        df[chrom_col] = df[chrom_col].astype(str)
        df[[hugo_col, patient_col, class_col]].head()
        if 'aa_change' in df.columns:
            df.rename(columns={'aa_change': 'aa_change_orig'}, inplace=True)
        df['aa_change'] = np.nan
        self.patient_col = patient_col
        self.sample_col = sample_col
        self.hugo_col = hugo_col
        self.chrom_col = chrom_col
        self.class_col = class_col
        self.ref_col = ref_col
        self.alt_col = alt_col
        self.start_col = start_col
        self.end_col = end_col
        return df

    def get_vep_df(self):
        new_names = OrderedDict([(self.chrom_col, 'Chromosome'),
                                 (self.start_col, 'Start_Position'),
                                 (self.ref_col, 'Reference_Allele'),
                                 (self.alt_col, 'Tumor_Seq_Allele2'),
                                 (self.sample_col, 'Tumor_Sample_Barcode'),
                                 ])

        m = self.df.loc[:, list(new_names.keys())].copy()
        m.rename(columns=new_names, inplace=True)
        m.Chromosome = m.Chromosome.apply(fix_chromosomes)
        return m


def get_patient_from_barcode(barcode):
    if barcode.startswith('TCGA'):
        return barcode[:12]
    else:
        return barcode


def fix_chromosomes(chrom):
    chrom = chrom.lstrip('chr')
    chrom = 'MT' if chrom == 'M' else chrom
    return chrom


if __name__ == "__main__" and __package__ is None:
    __package__ = "maf_matrix"
