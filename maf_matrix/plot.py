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

import matplotlib.collections as collections
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from collections import defaultdict

class MatrixPlotter(object):

    def __init__(self, matrix_df, n_all_patients=None):
        """Generate matrix plot using matrix dataframe.

        Args:
            matrix_df (pandas DataFrame): columns are hugo symbols, rows are
                patient_ids. Contents are boolean for mutation.
        """
        self.matrix_df = matrix_df
        self.use_patients = list(matrix_df.index)
        self.use_genes = list(matrix_df.columns)

        pbox_lgene = 25  # gene box width (pixels)
        pbox_lpatient = 25  # patient box width (pixels)
        self.dpi = 60

        n_genes = len(self.use_genes)
        n_patients = len(self.use_patients)

        # GENE LABELS
        if n_all_patients:
            p_counts = self.matrix_df.apply(sum, axis=0)
            max_chars = max([len(i) for i in self.use_genes])
            gene_labels = []
            format_str = "{gene:<" + str(max_chars) + "} {pc:3.0f}%"
            for gene in self.use_genes:
                gene_labels.append(format_str.format(
                    gene=gene, pc=float(p_counts[gene]) / n_all_patients * 100))
            n_mutated = sum(matrix_df.any(axis=1))
            pc_patients = n_mutated / float(n_all_patients) * 100
            self.title = "{n} patients ({pc:.1f}%)".format(n=n_mutated,
                                                           pc=pc_patients)
        else:
            gene_labels = self.use_genes
            self.title = "{n} patients".format(n=n_patients)

        self.ax_x = dict(length=pbox_lpatient * n_patients,
                         padding=90,
                         num=n_patients,
                         labels=self.use_patients,
                         boxlen=pbox_lpatient)
        self.ax_y = dict(length=pbox_lgene * n_genes,
                         padding=70,
                         num=n_genes,
                         labels=gene_labels,
                         boxlen=pbox_lgene)

    def draw_matrix(self):
        """Generate matrix figure."""
        figsize = [(self.ax_x['length'] + self.ax_x['padding']) / self.dpi,
                   (self.ax_y['length'] + self.ax_y['padding']) / self.dpi]
        fig1 = plt.figure(figsize=figsize)
        ax1 = fig1.add_subplot(111, aspect='equal')

        plt.xlim([0, self.ax_x['length']])
        plt.ylim([0, self.ax_y['length']])

        plt.gca().invert_yaxis()

        # ADD GRAY BACKGROUND
        ax1.add_patch(
            patches.Rectangle((0, 0), self.ax_x['length'], self.ax_y['length'],
                              fc=[0.85, 0.85, 0.85], ec='none')
        )
        # ADD MUTATION BOXES
        for p_ind, patient in enumerate(self.use_patients):
            for g_ind, gene in enumerate(self.use_genes):
                if self.matrix_df.iloc[p_ind, g_ind]:
                    self._add_box(p_ind, g_ind)
        # WHITE LINE GRID
        self._add_grid(ax1)

        ax1.set_frame_on(False)

        plt.xticks([self.ax_x['boxlen']*(d + 0.5)
                    for d in range(self.ax_x['num'])],
                   self.ax_x['labels'], ha='right',
                   rotation=45, rotation_mode="anchor")
        plt.yticks([self.ax_y['boxlen']*(d + 0.5)
                    for d in range(self.ax_y['num'])],
                   self.ax_y['labels'])

        ax1.tick_params(axis=u'both', which=u'both', length=0)
        ax1.annotate(self.title, xy=(1, 1), xycoords=ax1.transAxes,
                     ha='right', va='bottom')
        return fig1

    def _add_grid(self, ax):
        """Add white grid lines.

        One line 'after' each box excluding last one.

        Args:
            ax: Axis handle
            ax_x: dictionary of x-axis properties
            ax_y: dictionary of y-axis properties

        'labels', 'padding', 'length', 'num', 'boxlen'
        """
        n_x = self.ax_x['length'] / self.ax_x['boxlen']
        n_y = self.ax_y['length'] / self.ax_y['boxlen']
        # draw X grid
        if n_x > 1:
            xy_start = np.vstack((self.ax_x['boxlen'] *
                                  np.array(range(1, n_x, 1)),
                                  np.zeros(n_x - 1)))
            xy_end = np.vstack((self.ax_x['boxlen'] *
                                np.array(range(1, n_x, 1)),
                                self.ax_y['length'] * np.ones(n_x - 1)))
            segments = zip(zip(*xy_start), zip(*xy_end))
            lc_x = collections.LineCollection(segments, colors=(1, 1, 1, 1),
                                              linewidths=2)
            ax.add_collection(lc_x)
        # draw Y grid
        if n_y > 1:
            xy_start = np.vstack((np.zeros(n_y - 1),
                                  self.ax_y['boxlen'] *
                                  np.array(range(1, n_y, 1))))
            xy_end = np.vstack((self.ax_x['length'] * np.ones(n_y - 1),
                                self.ax_y['boxlen'] *
                                np.array(range(1, n_y, 1))))
            segments = zip(zip(*xy_start), zip(*xy_end))
            lc_y = collections.LineCollection(segments, colors=(1, 1, 1, 1),
                                              linewidths=2)
            ax.add_collection(lc_y)

        ax.add_patch(
            patches.Rectangle((0, 0), self.ax_x['length'], self.ax_y['length'],
                              fc='none', ec='white')
        )

    def _add_box(self, x_i, y_i):
        """Add rectangle for mutation using x-index and y-index."""
        l = x_i * self.ax_x['boxlen']
        b = y_i * self.ax_y['boxlen']
        w = self.ax_x['boxlen']
        h = self.ax_y['boxlen']
        rect = patches.Rectangle((l, b), w, h, fc='red', ec='none')
        plt.gca().add_patch(rect)
