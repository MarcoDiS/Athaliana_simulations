"""
19 Jul 2013
"""
from cPickle                          import load, dump
from subprocess                       import Popen, PIPE
from math                             import acos, degrees, pi, sqrt
from warnings                         import warn
from string                           import uppercase as uc, lowercase as lc
from random                           import random
from os.path                          import exists
from itertools                        import combinations
from uuid                             import uuid5, UUID
from hashlib                          import md5

from numpy                            import exp as np_exp
from numpy                            import median as np_median
from numpy                            import mean as np_mean
from numpy                            import std as np_std, log2
from numpy                            import array, cross, dot, ma, isnan
from numpy                            import histogram, linspace
from numpy                            import zeros
from numpy.linalg                     import norm

from scipy.optimize                   import curve_fit
from scipy.stats                      import spearmanr, pearsonr, chisquare
from scipy.stats                      import linregress
from scipy.stats                      import normaltest, norm as sc_norm
from scipy.cluster.hierarchy          import linkage, fcluster

from pytadbit                         import get_dependencies_version
from pytadbit.utils.three_dim_stats   import calc_consistency, mass_center
from pytadbit.utils.three_dim_stats   import dihedral, calc_eqv_rmsd
from pytadbit.utils.three_dim_stats   import get_center_of_mass, distance
from pytadbit.utils.tadmaths          import calinski_harabasz, nozero_log_list
from pytadbit.utils.tadmaths          import mean_none
from pytadbit.utils.extraviews        import plot_3d_model, setup_plot
from pytadbit.utils.extraviews        import chimera_view, tadbit_savefig
from pytadbit.utils.extraviews        import augmented_dendrogram, plot_hist_box
from pytadbit.utils.extraviews        import tad_coloring
from pytadbit.utils.extraviews        import tad_border_coloring
from pytadbit.utils.extraviews        import color_residues
from pytadbit.modelling.impmodel      import IMPmodel
from pytadbit.centroid                import centroid_wrapper
from pytadbit.aligner3d               import aligner3d_wrapper
from pytadbit.squared_distance_matrix import squared_distance_matrix_calculation_wrapper

import sys

try:
    from matplotlib import pyplot as plt
    from matplotlib.cm import jet, bwr
except ImportError:
    warn('matplotlib not found\n')

def contact_map(matrix, cutoff, nloci, axe=None):
    """
    Plots a contact map representing the frequency of interaction (defined
    by a distance cutoff) between two particles.
    
    :param None models: if None (default) the contact map will be computed
    using all the models. A list of numbers corresponding to a given set
    of models can be passed
    :param None cluster: compute the contact map only for the models in the
    cluster number 'cluster'
    :param None cutoff: distance cutoff (nm) to define whether two particles
    are in contact or not, default is 2 times resolution, times scale.
    :param None axe: a matplotlib.axes.Axes object to define the plot
    appearance
    :param None savefig: path to a file where to save the image generated;
    if None, the image will be shown using matplotlib GUI (the extension
    of the file name will determine the desired format).
    :param None savedata: path to a file where to save the contact map data
    generated, in three columns format (particle1, particle2, percentage
    of models where these two particles are in contact)
    
    """
    cmap = plt.get_cmap('jet')
    cmap.set_bad('darkgrey', 1)
    ims = axe.imshow(matrix, origin='lower', interpolation="nearest",
                     vmin=0, vmax=1, cmap=cmap,
                     extent=(0.5, nloci + 0.5, 0.5, nloci + 0.5))
    axe.set_ylabel('Particle')
    axe.set_xlabel('Particle')
    cbar = axe.figure.colorbar(ims)
    cbar.ax.set_yticklabels(['%3s%%' % (p) for p in range(0, 110, 10)])
    cbar.ax.set_ylabel('Percentage of models with particles at <' +
                       '%s nm' % (cutoff))
    axe.set_title('Contact map')


def correlate_with_real_data(model_contact_matrix, exp_interaction_matrix, nloci, cutoff, corr, 
                             off_diag, savefig, timestep, log_corr=True, axe=None):

    """
    Plots the result of a correlation between a given group of models and
    original Hi-C data.
    
    :param None models: if None (default) the correlation will be computed
    using all the models. A list of numbers corresponding to a given set
    of models can be passed
    :param None cluster: compute the correlation only for the models in the
    cluster number 'cluster'
    :param None cutoff: distance cutoff (nm) to define whether two particles
    are in contact or not, default is 2 times resolution, times scale.
    :param None savefig: path to a file where to save the image generated;
    if None, the image will be shown using matplotlib GUI (the extension
    of the file name will determine the desired format).
    :param False plot: to display the plot
    :param True log_corr: log plot for correlation
    :param None axe: a matplotlib.axes.Axes object to define the plot
    appearance
    :param None contact_matrix: input a contact matrix instead of computing
    it from the models
    
    :returns: correlation coefficient rho, between the two
    matrices. A rho value greater than 0.7 indicates a very good
    correlation
    
    """

    expdata = []
    moddata = []
    for i in xrange(len(exp_interaction_matrix)):
        for j in xrange(i + off_diag, len(exp_interaction_matrix)):
            #if not exp_interaction_matrix[i][j] > 0:
            #    continue
            expdata.append(exp_interaction_matrix[i][j])
            moddata.append(model_contact_matrix[i][j])

    if corr == 'spearman':
        corr = spearmanr(moddata, expdata)
    elif corr == 'pearson':
        corr = pearsonr(moddata, expdata)
    elif corr == 'logpearson':
        corr = pearsonr(nozero_log_list(moddata), nozero_log_list(expdata))
    elif corr == 'chi2':
        corr = chisquare(array(moddata), array(expdata))
        corr = 1. / corr[0], corr[1]
    else:
        raise NotImplementedError('ERROR: %s not implemented, must be one ' +
                                  'of spearman, pearson or frobenius\n')
    if not savefig:
        return corr
    if not axe:
        fig = plt.figure(figsize=(20, 4.5))
    else:
        fig = axe.get_figure()
    fig.suptitle('Correlation between normalized-real and modeled ' +
                 'contact maps (correlation=%.4f)' % (corr[0]),
                 size='x-large')
    ax = fig.add_subplot(131)

    # imshow of the modeled data
    #contact_map(model_contact_matrix, cutoff, nloci, axe=ax)
    cmap = plt.get_cmap('jet')
    cmap.set_bad('darkgrey', 1)
    ims = ax.imshow(log2(model_contact_matrix), origin='lower',
                    interpolation="nearest", cmap=cmap,
                    extent=(0.5, nloci + 0.5, 0.5, nloci + 0.5))
    ax.set_ylabel('Particles')
    ax.set_xlabel('Particles')
    ax.set_title('Models contacts count at timestep %s' % timestep)
    #ax.set_title('Model contacts count')
    cbar = ax.figure.colorbar(ims)
    cbar.ax.set_ylabel('Log2 (Model contacts count)')

    # correlation
    ax = fig.add_subplot(132)
    try:
        if log_corr:
            minmoddata = float(min([m for m in moddata if m]))
            minexpdata = float(min([m for m in expdata if m]))
            moddata, expdata = (log2([(m if m else minmoddata / 2) * 100 for m in moddata]),
                                log2([m if m else minexpdata / 2 for m in expdata]))
    except:
        warn('WARNING: unable to log for correlation with real data...')
    slope, intercept, r_value, p_value, _ = linregress(moddata, expdata)
    # slope, intercept, r_value, p_value, std_err = linregress(moddata, expdata)
    midplot = 'hexbin'
    if midplot == 'classic':
        lnr = ax.plot(moddata, intercept + slope * array (moddata), color='k',
                      ls='--', alpha=.7)
        ax.legend(lnr, ['p-value: %.3f, R: %.3f' % (p_value, r_value)])
        ax.plot(moddata, expdata, 'ro', alpha=0.5)
        ax.set_xlabel('Modelled data')
        ax.set_ylabel('Real data')
    elif midplot == 'hexbin':
        hb = ax.hexbin(moddata, expdata, mincnt=1,
                       gridsize=50, cmap=plt.cm.Spectral_r)
        lnr = ax.plot(moddata, intercept + slope * array (moddata), color='k',
                      ls='--', alpha=.7)
        ax.set_xlabel('Models contacts counts at dcutoff %d nm' % (cutoff))
        ax.set_ylabel('Normalized Hi-C count')
        cbaxes = fig.add_axes([0.41, 0.42, 0.005, 0.45])
        cbar = plt.colorbar(hb, cax=cbaxes)  # orientation='horizontal')
        cbar.set_label('Number of particle pairs')
    elif midplot == 'triple':
        maxval = max(expdata)
        minval = min(expdata)
        ax.set_visible(False)
        axleft = fig.add_axes([0.42, 0.18, 0.1, 0.65])
        axleft.spines['right'].set_color('none')
        axleft.spines['bottom'].set_color('none')
        axleft.spines['left'].set_smart_bounds(True)
        axleft.spines['top'].set_smart_bounds(True)
        axleft.xaxis.set_ticks_position('top')
        axleft.yaxis.set_ticks_position('left')
        axleft.set_ylabel('Normalized Hi-C count')
        axleft.patch.set_visible(False)
        axbott = fig.add_axes([0.44, 0.13, 0.17, 0.5])
        axbott.spines['left'].set_color('none')
        axbott.spines['top'].set_color('none')
        axbott.spines['left'].set_smart_bounds(True)
        axbott.spines['bottom'].set_smart_bounds(True)
        axbott.xaxis.set_ticks_position('bottom')
        axbott.yaxis.set_ticks_position('right')
        axbott.patch.set_visible(False)
        axbott.set_xlabel('Models contacts count')
        axmidl = fig.add_axes([0.44, 0.18, 0.17, 0.65])
        axbott.hist(moddata, bins=20, alpha=.2)
        x, _  = histogram([i if str(i) != '-inf' else 0. for i in expdata],
                          bins=20)
        axleft.barh(linspace(minval, maxval, 20), x,
                    height=(maxval - minval) / 20, alpha=.2)
        axleft.set_ylim((minval -
                         (maxval - minval) / 20, maxval +
                         (maxval - minval) / 20))
        axmidl.plot(moddata, expdata, 'k.', alpha=.3)
        axmidl.plot(moddata, intercept + slope * array (moddata), color='k',
                    ls='--', alpha=.7)
        axmidl.set_ylim(axleft.get_ylim())
        axmidl.set_xlim(axbott.get_xlim())
        axmidl.axis('off')
        # axmidl.patch.set_visible(False)
    ax.set_title('Real versus modelled data')
    ax = fig.add_subplot(133)
    cmap = plt.get_cmap('jet')
    cmap.set_bad('darkgrey', 1)
    ims = ax.imshow(log2(exp_interaction_matrix), origin='lower',
                    interpolation="nearest", cmap=cmap,
                    extent=(0.5, nloci + 0.5, 0.5, nloci + 0.5))
    ax.set_ylabel('Particles')
    ax.set_xlabel('Particles')
    ax.set_title('Normalized Hi-C count')
    cbar = ax.figure.colorbar(ims)
    cbar.ax.set_ylabel('Log2 (normalized Hi-C data)')
    if savefig:
        tadbit_savefig(savefig)
    elif not axe:
        plt.show()
    plt.close('all')
    return corr

# Getting input data
size     = int(sys.argv[1])
dcutoff  = float(sys.argv[2])
timestep = int(float(sys.argv[3])*0.006)

model_matrix = zeros((size,size))
exp_matrix   = zeros((size,size))

fpmodel = open("_model", "r")
for line in fpmodel.readlines():
    line  = line.strip()
    split = line.split()
    model_matrix[int(split[0])][int(split[1])] = float(split[2])
    model_matrix[int(split[1])][int(split[0])] = float(split[2])

fpexp = open("_exp", "r")
for line in fpexp.readlines():
    line  = line.strip()
    split = line.split()
    exp_matrix[int(split[0])][int(split[1])] = float(split[2])
    exp_matrix[int(split[1])][int(split[0])] = float(split[2])

print correlate_with_real_data(model_matrix, exp_matrix, nloci=size, cutoff=dcutoff, corr='spearman', off_diag=1, savefig="plot.png", timestep=timestep)[0]
