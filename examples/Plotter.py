#!/usr/bin/python3

import numpy as np
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
import pylab as plt



class Plotter():
    '''A generic (shell) class meant to streamline plotting. The intent is to modify this
    class for a specific plot.
    
    The structure of the class follows the convetions used in SciAnalysis:
    https://github.com/CFN-softbio/SciAnalysis/
    '''
    
    def __init__(self, x=None, y=None, name=None, plot_args=None, **kwargs):
        
        self.x = x
        self.y = y
        self.name = name

        self.x_label = kwargs['x_label'] if 'x_label' in kwargs else 'x'
        self.y_label = kwargs['y_label'] if 'y_label' in kwargs else 'y'
        self.x_rlabel = kwargs['x_rlabel'] if 'x_rlabel' in kwargs else self.x_label
        self.y_rlabel = kwargs['y_rlabel'] if 'y_rlabel' in kwargs else self.y_label
        self.x_err = kwargs['x_err'] if 'x_err' in kwargs else None
        self.y_err = kwargs['y_err'] if 'y_err' in kwargs else None
        
        self.plot_valid_keys = ['color', 'linestyle', 'linewidth', 'marker', 'markerfacecolor', 'markersize', 'alpha', 'markeredgewidth', 'markeredgecolor', 'capsize', 'ecolor', 'elinewidth']
        self.plot_args = { 'color' : 'k',
                        'marker' : 'o',
                        'linewidth' : 3.0,
                        'rcParams': {'axes.labelsize': 35,
                                        'xtick.labelsize': 30,
                                        'ytick.labelsize': 30,
                                        },
                            }        
        if plot_args: self.plot_args.update(plot_args)
        self._kwargs = kwargs # Save incase later methods depend on these settings
        
    
    def plot(self, save=None, show=False, plot_range=[None,None,None,None], plot_buffers=[0.2,0.05,0.2,0.05], **kwargs):
        '''Plots the data.
        
        Parameters
        ----------
        save : str
            Set to 'None' to avoid saving to disk. Provide filename to save.
        show : bool
            Set to true to open an interactive window.
        plot_range : [float, float, float, float]
            Set the range of the plotting (None scales automatically instead).
        '''  
        
        self._plot(save=save, show=show, plot_range=plot_range, plot_buffers=plot_buffers, **kwargs)
        
        
    def _plot(self, save=None, show=False, plot_range=[None,None,None,None], plot_buffers=[0.2,0.05,0.2,0.05], error=False, error_band=False, xlog=False, ylog=False, xticks=None, yticks=None, dashes=None, transparent=False, **kwargs):
        
        # DataLine._plot()
        
        plot_args = self.plot_args.copy()
        plot_args.update(kwargs)
        self.process_plot_args(**plot_args)
        
        
        
        self.fig = plt.figure( figsize=(10,7), facecolor='white' )
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )
        
        
        p_args = dict([(i, plot_args[i]) for i in self.plot_valid_keys if i in plot_args])
        self._plot_main(error=error, error_band=error_band, dashes=dashes, **p_args)
        
        
        self.ax.set_xlabel(self.x_rlabel)
        self.ax.set_ylabel(self.y_rlabel)
        
        if xlog:
            plt.semilogx()
        if ylog:
            plt.semilogy()
        if xticks is not None:
            self.ax.set_xticks(xticks)
        if yticks is not None:
            self.ax.set_yticks(yticks)

        if 'gridlines' in plot_args and plot_args['gridlines']:
            plt.grid()
        
        if 'title' in plot_args and plot_args['title'] is not None:
            size = plot_args['rcParams']['axes.labelsize']
            #size = plot_args['rcParams']['xtick.labelsize']
            plt.figtext(0, 1, plot_args['title'], size=size, weight='bold', verticalalignment='top', horizontalalignment='left')
        
        
        # Axis scaling
        xi, xf, yi, yf = self.ax.axis()
        if plot_range[0] != None: xi = plot_range[0]
        if plot_range[1] != None: xf = plot_range[1]
        if plot_range[2] != None: yi = plot_range[2]
        if plot_range[3] != None: yf = plot_range[3]
        self.ax.axis( [xi, xf, yi, yf] )
        
        self._plot_extra(**plot_args)
        
        if save:
            if 'dpi' in plot_args:
                plt.savefig(save, dpi=plot_args['dpi'], transparent=transparent)
            else:
                plt.savefig(save, transparent=transparent)
        
        if show:
            self._plot_interact()
            plt.show()
            
        plt.close(self.fig.number)
        
        

    def _plot_main(self, error=False, error_band=False, dashes=None, **plot_args):
        
        if error_band:
            # TODO: Make this work
            l, = plt.plot(self.x, self.y, **plot_args)
            self.ax.fill_between(self.x, self.y-self.y_err, self.y+self.y_err, facecolor='0.8', linewidth=0)
        
        elif error:
            l = plt.errorbar( self.x, self.y, xerr=self.x_err, yerr=self.y_err, **plot_args)
        
        else:
            #l, = plt.plot(self.x, self.y, **plot_args)
            l, = self.ax.plot(self.x, self.y, **plot_args)
            
            
        if dashes is not None:
            l.set_dashes(dashes)        
                
                
    def _plot_extra(self, **plot_args):
        '''This internal function can be over-ridden in order to force additional
        plotting behavior.'''
        
        pass
    
    

    def process_plot_args(self, **plot_args):
        
        if 'rcParams' in plot_args:
            for param, value in plot_args['rcParams'].items():
                plt.rcParams[param] = value    
                
                
                
if __name__ == '__main__':
    p = Plotter(x_rlabel='$x \, (\mu \mathrm{m})$', y_rlabel='$T \, (\mathrm{^{\circ}C})$')
    p.x = np.linspace(0, 100, num=100)
    p.y = np.square(p.x)

    p.plot(save='output.png')
