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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.collections as collections
from matplotlib.transforms import Bbox


__author__ = 'Stephen G. Gaffney'

COLOR_DICT = {
    'missense_mutation': (0.0, 0.4470588235294118, 0.6980392156862745),
    'nonsense_mutation': (0.0, 0.6196078431372549, 0.45098039215686275),
    'nonstop_mutation':
        (0.33725490196078434, 0.7058823529411765, 0.9137254901960784),
    'rna': (0.33725490196078434, 0.7058823529411765, 0.9137254901960784),
    'silent': (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
    'splice_site': (0.8352941176470589, 0.3686274509803922, 0.0),
    'translation_start_site':
        (0.33725490196078434, 0.7058823529411765, 0.9137254901960784),
    'other':
        (0.8, 0.4745098039215686, 0.6549019607843137)
}


class MatrixPlotter(object):
    def __init__(self, matrix_df, info_df=None, n_all_patients=None,
                 hypermutated=None, color_dict=None, show_title=True,
                 show_all_patients=False):
        """Generate matrix plot using matrix dataframe.

        Args:
            matrix_df (pandas DataFrame): columns are hugo symbols, rows are
                patient_ids. Contents are boolean for mutation, or string for
                annotation.
            info_df (pandas DataFrame). optional.
        """
        if not show_all_patients:
            n = matrix_df.apply(lambda s: (s != False).any()).sum()
            matrix_df = matrix_df.iloc[:, :n]

        self.matrix_df = matrix_df
        self.info_df = info_df
        self.color_dict = color_dict if color_dict else dict()
        self.show_title = show_title

        # self.m_status = ~(matrix_df == False)
        self.m_status = matrix_df.applymap(lambda x: x is not False)
        # assert isinstance(self.m_status, pd.DataFrame)
        self.use_genes = list(matrix_df.index)
        self.use_patients = list(matrix_df.columns)
        self.hypermutated = hypermutated if hypermutated is not None else []

        self.upi = float(96)  # units per inch

        self.n_genes = len(self.use_genes)
        self.n_patients = len(self.use_patients)

        # GENE LABELS
        if n_all_patients:
            p_counts = self.m_status.sum(axis=1)
            gene_labels = []
            format_str = "{gene} {pc:3.0f}%"
            for gene in self.use_genes:
                gene_labels.append(format_str.format(
                    gene=gene, pc=float(p_counts[gene]) / n_all_patients * 100))
            self.gene_labels = gene_labels
            n_mutated = sum(matrix_df.any(axis=0))
            pc_patients = n_mutated / float(n_all_patients) * 100
            patient_str = 'patient' if self.n_patients == 1 else 'patients'
            self.title = "{n} {patients} ({pc:.1f}%)".format(n=n_mutated,
                                                             patients=patient_str,
                                                             pc=pc_patients)
        else:
            self.gene_labels = self.use_genes
            self.title = "{n} patients".format(n=self.n_patients)
        patient_labels = []
        for i in self.use_patients:
            label = '{}*'.format(i) if i in self.hypermutated else i
            patient_labels.append(label)
        self.patient_labels = patient_labels

    def draw_matrix(self, fontsize=10, max_label=100, box_px=30,
                    show_limits=False):
        """Generate matrix figure."""
        # figsize = [(self.ax_x['length'] + self.ax_x['padding']) / self.upi,
        #            (self.ax_y['length'] + self.ax_y['padding']) / self.upi]

        pbox_lgene = float(box_px)  # gene box width (data units)
        pbox_lpatient = float(box_px)  # patient box width (data units)

        self.ax_x = dict(length=pbox_lpatient * self.n_patients,
                         # padding=float(90),
                         num=self.n_patients,
                         labels=self.patient_labels,
                         boxlen=pbox_lpatient)
        self.ax_y = dict(length=pbox_lgene * self.n_genes,
                         # padding=float(70),
                         num=self.n_genes,
                         labels=self.gene_labels,
                         boxlen=pbox_lgene)

        figsize = [self.ax_x['length'] / self.upi,
                   self.ax_y['length'] / self.upi]
        subplot_kw = {'aspect': 'equal', 'position': [0, 0, 1, 1],
                      'xlim': [0, self.ax_x['length']],
                      'ylim': [0, self.ax_y['length']]}
        hfig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw=subplot_kw)
        ax.invert_yaxis()
        obj_list = []  # for calculating figure extents
        # ADD GRAY BACKGROUND (orig [0.85, 0.85, 0.85])
        rect = patches.Rectangle((0, 0), self.ax_x['length'],
                                 self.ax_y['length'], fc='#EAEAF2',
                                 ec='none')
        patch = ax.add_patch(rect)
        obj_list.append(patch)
        # ADD MUTATION BOXES
        m_reindex = self.m_status.set_index(np.arange(0, self.n_genes))
        m_reindex.columns = np.arange(0, self.n_patients)
        s = m_reindex.stack()
        mut_inds = s[s].index.values
        for g_ind, p_ind in mut_inds:
            gene = self.use_genes[g_ind]
            patient = self.use_patients[p_ind]
            if self.info_df is not None:
                info_dict = self.info_df.loc[(patient, gene)].to_dict()
                if self.color_dict:
                    try:
                        color = self.color_dict[info_dict['var_class']]
                    except KeyError:
                        color = COLOR_DICT['other']
                else:
                    color = 'r'
                annot = info_dict['annot']
            else:
                color = 'r'
                annot = None

            self._add_box(p_ind, g_ind, ax, color=color)
            # annot = self.matrix_df.iloc[g_ind, p_ind]
            if type(annot) == str:
                self._add_annot(p_ind, g_ind, ax, annot, fontsize=fontsize,
                                fontweight='bold')
        # WHITE LINE GRID
        self._add_grid(ax)

        ax.set_frame_on(False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.tick_params(axis=u'both', which=u'both', length=0)

        labels = self._add_ax_labels(fontsize=fontsize)
        obj_list.extend(labels)

        # ADJUST FIGURE
        bb = get_target_bbox(obj_list, ax, hfig)
        # print(bb)

        # RIGHT WHITESPACE PADDING FOR KNOWN WIDTH
        if max_label is not None:
            title_overflow = max(bb.xmax - self.ax_x['length'], 0)
            r_padding = max([max_label + bb.xmin - title_overflow, 0])
        else:
            r_padding = 0
        # VERTICAL PADDING TO AVOID CROPPING TEXT
        v_padding = 5
        h_padding = 5

        # print("R padding: {}".format(r_padding))
        xmin = bb.xmin - h_padding
        xmax = bb.xmax + h_padding + r_padding
        ymin = bb.ymin - v_padding
        ymax = bb.ymax + v_padding

        bb_new = Bbox(np.array([[xmin, ymax], [xmax, ymin]]))
        # print(bb_new)
        if show_limits:
            rect = patches.Rectangle(bb_new.p0, bb_new.width, bb_new.height,
                                     alpha=0.2)
            ax.add_artist(rect)

        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymax, ymin])

        fwidth_in = (xmax - xmin) / self.upi
        fheight_in = (ymax - ymin) / self.upi
        hfig.set_figwidth(fwidth_in)
        hfig.set_figheight(fheight_in)

        return hfig, ax

    def _add_grid(self, ax):
        """Add white grid lines.

        One line 'after' each box excluding last one.
        ax_x/y has 'labels', 'padding', 'length', 'num', 'boxlen'

        Args:
            ax: Axis handle
        """
        n_x = int(self.ax_x['length'] / self.ax_x['boxlen'])
        n_y = int(self.ax_y['length'] / self.ax_y['boxlen'])
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

    def _add_box(self, x_i, y_i, ax, color='red'):
        """Add rectangle for mutation using x-index and y-index."""
        l = x_i * self.ax_x['boxlen']
        b = y_i * self.ax_y['boxlen']
        w = self.ax_x['boxlen']
        h = self.ax_y['boxlen']
        rect = patches.Rectangle((l, b), w, h, fc=color, ec='none')
        ax.add_patch(rect)

    def _add_annot(self, x_i, y_i, ax, annot, fontsize=8, fontweight='bold',
                   color='w'):
        """Add annotation string to mutation box using x,y index."""
        x = x_i * self.ax_x['boxlen'] + self.ax_x['boxlen'] / 2
        y = y_i * self.ax_y['boxlen'] + self.ax_y['boxlen'] / 2
        t = ax.text(x, y, annot, ha='center', va='center',
                    fontsize=fontsize, color=color, fontweight=fontweight)
        return t

    def _add_ax_labels(self, fontsize=10):
        """Add manual x and y labels"""
        labels = []
        for ind, ylabel in enumerate(self.ax_y['labels']):
            x = 0
            y = ind * self.ax_y['boxlen'] + self.ax_y['boxlen'] / 2
            t = plt.text(x, y, ylabel + ' ', ha='right', va='center',
                         fontsize=fontsize)
            labels.append(t)
        for ind, xlabel in enumerate(self.ax_x['labels']):
            x = ind * self.ax_x['boxlen'] + self.ax_x['boxlen'] / 2
            y = self.ax_y['length']
            t = plt.text(x, y, xlabel + ' ', ha='center', va='top',
                         fontsize=fontsize, rotation=90,
                         rotation_mode="default"
                         )
            labels.append(t)
        # Add title
        x = 0  # self.ax_x['length']
        y = 0
        if self.show_title:
            t = plt.text(x, y, self.title, ha='left', va='bottom',
                         fontsize=fontsize)
            labels.append(t)
        return labels


def get_target_bbox(obj_list, ax=None, hfig=None):
    """Get bounding box for iterable of plot objects.

    SEE:
    # BBOX: https://stackoverflow.com/questions/24581194/matplotlib-text-bounding-box-dimensions
    # https://stackoverflow.com/questions/22667224/matplotlib-get-text-bounding-box-independent-of-backend
    """
    transf = ax.transData.inverted()
    renderer = hfig.canvas.get_renderer()
    bb_list = []
    for obj in obj_list:
        bb_list.append(
            obj.get_window_extent(renderer=renderer).transformed(transf))
    extents = np.vstack([bb.extents for bb in bb_list])
    # print(extents)

    x0 = max(extents[:, 0]) if ax.xaxis_inverted() else min(extents[:, 0])
    y0 = max(extents[:, 1]) if ax.yaxis_inverted() else min(extents[:, 1])
    x1 = min(extents[:, 2]) if ax.xaxis_inverted() else max(extents[:, 2])
    y1 = min(extents[:, 3]) if ax.yaxis_inverted() else max(extents[:, 3])

    bbox = Bbox(np.array([[x0, y0], [x1, y1]]))
    # print(bbox)
    return bbox
