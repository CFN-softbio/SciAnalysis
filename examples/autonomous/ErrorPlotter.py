#!/usr/bin/python3

import numpy as np
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
import pylab as plt

# Helpers
########################################
# TODO: DEPRECATED in favor of tools.val_stats
def print_d(d, i=4):
    '''Simple helper to print a dictionary.'''
    for k, v in d.items():
        if isinstance(v,dict):
            print('{}{} : <dict>'.format(' '*i,k))
            print_d(v, i=i+4)
        elif isinstance(v,(np.ndarray)):
            print('{}{} : Ar{}: {}'.format(' '*i,k,v.shape,v))
        elif isinstance(v,(list,tuple)):
            print('{}{} : L{}: {}'.format(' '*i,k,len(v),v))
        else:
            print('{}{} : {}'.format(' '*i,k,v))

def print_results(results):
    '''Simple helper to print out a list of dictionaries.'''
    for i, result in enumerate(results):
        print(i)
        print_d(result)

def print_n(d):
    '''Simple helper to print nested arrays/dicts'''
    if isinstance(d, (list,tuple,np.ndarray)):
        print_results(d)
    elif isinstance(d, dict):
        print_d(d)
    else:
        print(d)

def val_stats(values, name='z'):
    span = np.max(values)-np.min(values)
    print("  {} = {:.2g} ± {:.2g} (span {:.2g}, from {:.3g} to {:.3g})".format(name, np.average(values), np.std(values), span, np.min(values), np.max(values)))



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
                        'marker' : None,
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
            size = plot_args['rcParams']['axes.labelsize']*0.5
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
        
        #l, = plt.plot(self.x, self.y, **plot_args)
        l, = self.ax.plot(self.x, self.y, alpha=0.25, **plot_args)
        
        from scipy.ndimage import gaussian_filter1d
        
        y_smoothed = gaussian_filter1d(self.y, sigma=10, mode='nearest')
        
        l, = self.ax.plot(self.x, y_smoothed, **plot_args)
        if dashes is not None:
            l.set_dashes(dashes)        


        N = int(len(self.y)*0.10)
        y_final = np.average(self.y[-N:])
        self.ax.axhline(y_final, color='b', alpha=0.5, linewidth=1.0)
        
        thresh = 0.8
        y_span = self.y[0]-y_final
        y_target = self.y[0] - y_span*thresh
        xt, yt = self.target_y(np.asarray(self.x), np.asarray(y_smoothed), y_target)
        self.ax.plot(xt, yt, 'o', color='b', markersize=10, alpha=0.75)
        s = '$\epsilon$ decrease reaches {:.0f}% by $N = {:d}$'.format(100*thresh, xt)
        e = y_span*0.3
        self.ax.text(xt, yt+e, s, size=15, color='b', verticalalignment='bottom', horizontalalignment='left')
        self.ax.plot([xt, xt], [yt, yt+e], color='b', alpha=0.5, linewidth=0.5)
        
        
    def target_y(self, x, y, target):
        '''Find the datapoint closest to the given y.'''
    
        #x = np.asarray(self.x)
        #y = np.asarray(self.y)

        # Sort
        indices = np.argsort(y)
        x_sorted = x[indices]
        y_sorted = y[indices]

        # Search through y for the target
        idx = np.where( y_sorted>=target )[0][0]
        xcur = x_sorted[idx]
        ycur = y_sorted[idx]
        
        return xcur, ycur
        
            
                
                
    def _plot_extra(self, **plot_args):
        '''This internal function can be over-ridden in order to force additional
        plotting behavior.'''
        
        pass
    
    

    def process_plot_args(self, **plot_args):
        
        if 'rcParams' in plot_args:
            for param, value in plot_args['rcParams'].items():
                plt.rcParams[param] = value    
                
                
                
if __name__ == '__main__':
    
    
    #infile = 'Data_2021-02-09_10_33_48_model_1.npy' # 101
    #infile, sample_name = 'Data_2021-02-09_11_59_02_model_1.npy', 'sample2_anneal15_run1' # 1341
    #infile, sample_name = 'Data_2021-02-09_21_54_05_model_1.npy', 'sample3_anneal1200_run1' # 1502
    #infile, sample_name = 'Data_2021-02-10_08_57_09_model_1.npy', 'sample5_anneal300_run1' # 1312 
    infile, sample_name = 'Data_2021-02-10_18_29_11_model_1.npy', 'sample1_anneal5_run1' # 1244 
    #infile = 'Data_2021-02-11_03_39_16_model_1.npy' # 4
    #infile = 'Data_2021-02-11_03_43_42_model_1.npy' # 767
    
    data = np.load(infile, allow_pickle=True)
    print('Loaded {} records from {}'.format(data.shape[0], infile))
    print('    filename: {}'.format(data[-1]['filename']))
    
    vals = np.asarray([d['time stamp'] for d in data])
    idx = np.argsort(vals)
    data = data[idx]
    
    #vals = np.asarray([d['objective function evaluation'] for d in data])
    vals = []
    for d in data:
        if 'objective function evaluation' in d:
            vals.append(d['objective function evaluation'])
    Ns = range(len(vals))
    
    
    p = Plotter(x=Ns, y=vals, x_rlabel='$N$', y_rlabel='$ \epsilon \, (\mathrm{a.u.})$')

    p.plot(save='error_history-{}.png'.format(sample_name), title=sample_name, plot_range=[0, None, None, None], plot_buffers=[0.14, 0.04, 0.15, 0.08])
