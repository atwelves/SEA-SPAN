# Modified further by Andrew G Twelves
# Last edit 21.10.26

# -*- coding: utf-8 -*-
#
#****************************************************************************#
#                                                                            #
#  Copyright (c) 2015, by University of Birmingham. All rights reserved.     #
#                                                                            #
#  Redistribution and use in source and binary forms, with or without        #
#  modification, are permitted provided that the following conditions        #
#  are met:                                                                  #
#                                                                            #
#      * Redistributions of source code must retain the above copyright      #
#        notice, this list of conditions and the following disclaimer.       #
#      * Redistributions in binary form must reproduce the above copyright   #
#        notice, this list of conditions and the following disclaimer in the #
#        documentation and/or other materials provided with the distribution.#
#      * The name 'University of Birmingham' may not be used to endorse or   #
#        promote produces derived from this software without specific prior  #
#        written permission.                                                 #
#                                                                            #
#****************************************************************************#
#
"""
Modified Taylor diagram [Elvidge 2014, Taylor 2001] implementation code.

A set of classes and functions for constructiong a "modified Taylor diagram".
A diagram which plots the model standard deviation, correlation,
bias, error standard deviation and mean square error to observation data
with r=stddev and theta=arccos(correlation).

This code can be used for creating your own diagrams, and an example use
is shown at the bottom of the program. If you use this approach for model
comparison we ask that you reference:

Elvidge, S., M. J. Angling, and B. Nava (2014), On the Use of Modified Taylor 
Diagrams to Compare Ionospheric Assimilation Models, Radio Sci., 
doi:10.1002/2014RS005435.

The orignal Taylor paper is:
Taylor, K. E. (2001), Summarizing multiple aspects of model performance in a 
single diagram, J. Geophys. Res., 106(7), 7183–7192.

The format of use is (a more detailed example can be found at bottom of 
the file):

>>> dia = ModifiedTaylorDiagram(refstd, label="Observation")
>>> colors = dia.calc_colors(bias_data)
>>> cbar = dia.add_colorbar(fig, bias_data.min(), bias_data.max())
>>> contours = dia.add_contours(colors='darkblue', levels=np.linspace(0,1.5,7))
>>> plt.clabel(contours, contours.levels[0:len(contours.levels)-2],
               inline=1, fontsize=10, use_clabeltext=1)
>>> plt.text(0.1, 0.9, 'Error Std. Dev.', transform=dia.ax.transAxes, 
             rotation=30, color='darkblue', fontsize=10)
>>> for i,(stddev,corrcoef,bias) in enumerate(data):
        dia.add_point(stddev, corrcoef, mod_lab[i], marker='o', markersize=10., 
                      ls='', c=cm.jet(colors[i]), label=labels)
>>> legend = dia.add_legend(str_artist, labels, prop=dict(size='small'), 
                            loc=(0.8,0.9))
>>> plt.show()

The code is split up in such as way to give the user as much flexible in
presentation (color etc.) as wanted.

Modification History
-------
Created on Mon Nov 01 2014 by Sean at SERENE, University of Birmingham
Contact: spaceweather@contacts.bham.ac.uk

This code is an expansion of an original piece of work by Yannick Copin 
(yannick.copin@laposte.net) created on 2012-02-17.

03/08/15  Added 'offset' variable to add_point, and 'fontsize' to add_colorbar
10/12/24  Completely updated code
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.text import Text
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf
from matplotlib.projections import PolarAxes
from matplotlib.ticker import MaxNLocator
import cmocean

import xarray as xr

class ArtistObject:
    def __init__(self, text):
        self.my_text = text


class LegendHandler:
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        # The variables x0, y0, width, height were not needed
        patch = Text(x=0, y=0, text=orig_handle.my_text,
                     verticalalignment='baseline',
                     horizontalalignment='left')
        handlebox.add_artist(patch)
        return patch


class ModifiedTaylorDiagram:
    def __init__(self, datastd, fig=None, rect=111, normalize='Y',
                 sd_axis_frac=1.5, fontsize=10):

        # Normalization logic
        if normalize.upper() == 'Y':
            self.datastd = 1.0
            self.normfactor = datastd
        else:
            self.datastd = datastd
            self.normfactor = 1.0

        # Save figure if passed
        if fig is None:
            fig = plt.figure()
        self._fig = fig

        self.smin = 0
        self.smax = sd_axis_frac * self.datastd

        tr = PolarAxes.PolarTransform()

        # Correlation locations and tick formatting
        rlocs = np.concatenate((np.arange(10) / 10., [0.95, 0.99]))
        tlocs = np.arccos(rlocs)
        gl1 = gf.FixedLocator(tlocs)
        tf1 = gf.DictFormatter({tl: f"{rl}" for tl, rl in zip(tlocs, rlocs)})

        # Std. Dev. locations and formatting
        # Using a stable approach for rounding
        sdlocs_all = np.linspace(self.smin, self.smax, 7)
        # Round values nicely
        sdlocs = [float(f"{val:.2g}") for val in sdlocs_all[1:]]
        sdlocs.append(0)
        sdlocs = sorted(sdlocs)  # Ensure ascending order if needed
        gl2 = gf.FixedLocator(sdlocs)
        tf2 = gf.DictFormatter({sd: f"{sd}" for sd in sdlocs})

        ghelper = fa.GridHelperCurveLinear(tr,
                                           extremes=(0, np.pi / 2, self.smin, self.smax),
                                           grid_locator1=gl1,
                                           grid_locator2=gl2,
                                           tick_formatter1=tf1,
                                           tick_formatter2=tf2
                                           )

        ax = fa.FloatingSubplot(self._fig, rect, grid_helper=ghelper)
        self._fig.add_subplot(ax)

        # Set axis directions
        ax.axis["top"].set_axis_direction("bottom")
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].major_ticklabels.set_color("k")
        ax.axis["top"].label.set_color("k")
        ax.axis["top"].label.set_text("Correlation")
        ax.axis["top"].label.set_fontsize(20)
        ax.axis["top"].major_ticklabels.set_fontsize(12)

        ax.axis["left"].set_axis_direction("bottom")
        ax.axis["left"].toggle(ticklabels=True)
        ax.axis["left"].major_ticklabels.set_fontsize(12)

        ax.axis["right"].set_axis_direction("top")
        ax.axis["right"].toggle(ticklabels=True, label=True)
        ax.axis["right"].major_ticklabels.set_axis_direction("left")
        ax.axis["right"].major_ticklabels.set_fontsize(12)
        if normalize.upper() == 'Y':
            ax.axis["right"].label.set_text("Normalised standard deviation")
        else:
            ax.axis["right"].label.set_text("Standard deviation")
        ax.axis["right"].label.set_fontsize(20)

        ax.axis["bottom"].set_visible(False)

        ax.grid(False)

        self._ax = ax  # The floating subplot
        self.ax = ax.get_aux_axes(tr)  # Polar coordinates

        # Reference point and std. dev. contour
        l, = self.ax.plot([0], self.datastd, marker=".", markerfacecolor='k', markeredgecolor='k', ls='', ms=12)
        t = np.linspace(0, np.pi / 2)
        r = np.zeros_like(t) + self.datastd
        self.ax.plot(t, r, 'k--', label='_', alpha=0.75)

        self.samplePoints = [l]

        # Add correlation lines for reference
        for i in np.concatenate((np.arange(1, 10) / 10., [0.99])):
            self.ax.plot([np.arccos(i), np.arccos(i)], [0, self.smax],
                         c='k', alpha=0.35)

    def add_point(self, obsstd, modstd, corrcoef, modlab, offset=(-12, -12), **kwargs):
        print(type(modstd))
        print(type(obsstd))
        modstd = modstd / obsstd
        l, = self.ax.plot(np.arccos(corrcoef), modstd, **kwargs)
        self.samplePoints.append(l)
        #self.ax.annotate(modlab, xy=(np.arccos(corrcoef), modstd),
        #                 xytext=offset, textcoords='offset points')
        return l

    def add_contours(self, levels=6, **kwargs):
        rs, ts = np.meshgrid(np.linspace(self.smin, self.smax),
                             np.linspace(0, np.pi / 2))
        rms = np.sqrt(self.datastd ** 2 + rs ** 2 - 2 * self.datastd * rs * np.cos(ts))
        contours = self.ax.contour(ts, rs, rms, levels, alpha=0.75, **kwargs)
        return contours

    def calc_colors(self, bias):
        bias = bias / self.normfactor
        largest = max(abs(bias.min()), abs(bias.max()))
        bias_col = (bias + largest) / (2 * largest)
        return bias_col

    def add_colorbar(self, bias, num_ticks=5, fontsize=10):
        bias = bias / self.normfactor
        cnorm = Normalize(vmin=0, vmax=15)

        sm = plt.cm.ScalarMappable(norm=cnorm, cmap=cm.cool)
        sm.set_array([])

        # Explicitly reference the figure and the main axes for colorbar
        cbar = self._fig.colorbar(sm, ax=self._ax, orientation='horizontal',
                                  fraction=0.042, pad=0.1)
        cbar.locator = MaxNLocator(nbins=num_ticks)
        cbar.update_ticks()
        cbar.outline.set_visible(False)
        cbar.set_label('Scale', fontsize=fontsize)

    def add_legend(self, artist, labels, *args, **kwargs):
        obj = [ArtistObject(st) for st in artist]
        if self.normfactor != 1:
            if self.normfactor > 0.01:
                obj.append(ArtistObject(f"{self.normfactor:4.2f}"))
                labels.append('Norm factor')
            else:
                obj.append(ArtistObject(f"{self.normfactor:4.2e}"))
                labels.append('Norm factor')

        dic = {o: LegendHandler() for o in obj}
        self.ax.legend(obj, labels, handler_map=dic, *args, **kwargs)


if __name__ == '__main__':
    
    stations = ['Kemi','Oulu','Raahe',
            'Pietarsaari','Vaasa',
            'Kaskinen','Pori','Rauma',
            'Degerby',
            'Hanko','Helsinki','Porvoo','Hamina']

#    stations = ['Degerby']

    # set colormap - different colours for Gulf of Bothnia, Gulf of Finland and Degerby
    custom = [( 27/255,158/255,119/255),( 27/255,158/255,119/255),( 27/255,158/255,119/255),
              ( 27/255,158/255,119/255),( 27/255,158/255,119/255),
              ( 27/255,158/255,119/255),( 27/255,158/255,119/255),( 27/255,158/255,119/255),
              (217/255, 95/255,  2/255),
              (117/255,112/255,179/255),(117/255,112/255,179/255),(117/255,112/255,179/255),(117/255,112/255,179/250)]

    markers = ['v','^','>',
           '<','8',
           'p','s','h',
           'o',
           '*','P','X','D'] 

    fig = plt.figure(figsize=(11, 7))
    dia = ModifiedTaylorDiagram(1, fig=fig, normalize='Y', fontsize=14)

    # Choose station
    st = 'Kemi'

    # Declare empty arrays
    obs = np.zeros((16,366))
    mod = np.zeros((16,366))

    # Open mesoscale part of observed spectrum
    ds  = xr.open_mfdataset('{}_obs_meso_*'.format(st))    
    mes = ds.__xarray_dataarray_variable__.values
    obs[0:8 ,:] = np.transpose(mes[:,:])
    
    # Open synoptic-scale part of observed spectrum
    ds  = xr.open_mfdataset('{}_obs_syno_*'.format(st))
    syn = ds.__xarray_dataarray_variable__.values
    obs[8:16,:] = np.transpose(syn[:,:])

    # Open mesoscale part of model spectrum
    ds = xr.open_mfdataset('{}_mod_meso_*'.format(st))
    mes = ds.__xarray_dataarray_variable__.values
    mod[0:8,:] = np.transpose(mes[:,:])

    # Open synoptic-scale part of model spectrum
    ds = xr.open_mfdataset('{}_mod_syno_*'.format(st))
    syn = ds.__xarray_dataarray_variable__.values
    mod[8:16,:] = np.transpose(syn[:,:])

    # Calculate sum of bias for scaling of scatter points
    sum_bias = np.nansum(mod) - np.nansum(obs)
    print(sum_bias)

    for f in range(0,15):
        print(f)
        fstr = str(f)
        data = obs[f,:]
        m    = mod[f,:]

        # Calculate standard deviation
        datastd = data.std(ddof=1)
        points = np.array([m.std(ddof=1), np.corrcoef(data, m)[0, 1], m.mean() - data.mean()])        

        # Standard deviation in spectral power
        modstd = points[0]
        # Correlation in spectral power
        corrcoef = points[1]
        # Bias in spectral power
        bias   = 366*points[2]

        print(bias)
        dia.add_point(datastd,modstd, corrcoef,fstr, marker='o',
                      markersize=300*bias/sum_bias,c=cm.cool(f/15),alpha = 0.5)

    dia.add_legend(stations, [f"Model {lbl}" for lbl in stations],
                   prop=dict(size='small'), loc=(0.8, 0.9))

    contours = dia.add_contours(colors='k',
                                levels=np.linspace(0, 1.5, 7),
                                linestyles='dotted')
    plt.clabel(contours, inline=1, use_clabeltext=1, fmt='%1.2f',fontsize=12)
    dia.ax.text(1.5, 1.25, 'Error Std. Dev.',rotation=30, color='k',fontsize=12)

    fig.savefig("fig.png", dpi=600, facecolor='w', edgecolor='w')
    plt.savefig('taylor_wavelet_{}.png'.format(st))

    plt.figure()
    plt.pcolormesh([[2,9.5],[2,9.5]],cmap=cm.cool)
    cbar = plt.colorbar()
    cbar.set_ticks([2,4,6,8],labels=['4','16','64','256'],fontsize=20)
    plt.savefig('taylor_comp_legend.png')
