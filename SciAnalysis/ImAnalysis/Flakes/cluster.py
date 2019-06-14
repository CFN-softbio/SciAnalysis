#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4

import pickle
from ..Protocols import *


from scipy.spatial.distance import cdist
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, MeanShift, estimate_bandwidth, AffinityPropagation # Clustering methods
from sklearn.decomposition import PCA

from mpl_toolkits.mplot3d import Axes3D

class cluster(ProtocolMultiple):
    
    def __init__(self, name='cluster', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.pkl'
        self.run_args = {
                        'file_extension' : '.pkl',
                        'force' : False,
                        'verbosity' : 3,
                        'num_jobs' : None,
                        'num_clusters' : 20,
                        'cluster_method' : 'kmeans', # 'affinity', 'kmeans', 'meanshift'
                        'features' : 'all', # 'shape', 'color', 'all'
                        }
        self.run_args.update(kwargs)
        
        
    def load_flakes(self, datas, **run_args):
        flakes = []
        for data in datas:
            with open(data.infile, 'rb') as fin:
                saved = pickle.load(fin) # 'res_map', 'image_labelmap', 'flakes'
                if len(flakes)==0:
                    h, w = saved['res_map'].shape
                
                for flake in saved['flakes']:
                    flakes.append(flake)
                    #self.print_structure(flake)
                    
                if run_args['verbosity']>=5:
                    print('      {} flakes added from image {}'.format(len(saved['flakes']), data.infile))        
                    
        return flakes
    
    
    def load_flakes_parallel(self, datas, **run_args):
        # Parallelize loading
        # Doesn't seem to actually run faster (likely I/O limited)
        from joblib import Parallel, delayed
        import itertools
        
        flakes = Parallel(n_jobs=run_args['num_jobs'])(delayed(self.load_flake_pkl)(data.infile) for data in datas)
        flakes = list(itertools.chain.from_iterable(flakes))
        
        return flakes
    
    def load_flake_pkl(self, infile):
        with open(infile, 'rb') as fin:
            flakes = pickle.load(fin)['flakes']
        return flakes
    

    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}
        to_save = {}
        
        # Aggregate results
        ########################################
        flakes = self.load_flakes(datas, **run_args)
        if run_args['verbosity']>=4:
            print('  {:,d} flakes identified in {:d} images'.format(len(flakes), len(datas)))
            
        

        if run_args['features']=='all':
            features = [ np.concatenate([flake['flake_color_fea'], flake['flake_shape_fea']]) for flake in flakes ]
            
        else:
            features = [ flake['flake_{}_fea'.format(run_args['features'])] for flake in flakes ]
        
        features_orig = np.asarray(features)
        
        
        
        # Clustering
        ########################################
        features = StandardScaler().fit_transform(features_orig)
        
        if run_args['verbosity']>=4:
            print("  Clustering {:,d} flakes using '{}'".format(len(flakes), run_args['cluster_method']))
            
        start = time.time()
                                                             
        if run_args['cluster_method']=='kmeans':
            cluster_result = KMeans(n_clusters=run_args['num_clusters'], random_state=0, n_jobs=-1).fit(features)
            
        elif run_args['cluster_method']=='meanshift':
            bandwidth = estimate_bandwidth(features, quantile=0.1)#, n_samples=int(features.shape[0]/10))
            cluster_result = MeanShift(bandwidth=bandwidth, bin_seeding=True).fit(features)
            
        elif run_args['cluster_method']=='affinity':
            cluster_result = AffinityPropagation().fit(features)
            
        else:
            print("ERROR: clustering method '{}' not recognized.".format(run_args['cluster_method']))
            raise NotImplementedError

        to_save['cluster_result'] = cluster_result
        results['cluster_runtime'] = time.time()-start
        results['cluster_method'] = run_args['cluster_method']
            
        # Assignments are unsorted by default
        assignment = cluster_result.labels_
        results['num_clusters'] = len(np.unique(assignment))

        if run_args['verbosity']>=4:
            print("    clustering took {:.1f}s ({:d} clusters)".format(results['cluster_runtime'], results['num_clusters']))
            
            
        # Sort clusters into a sensible order
        consider_features = np.asarray([flake['flake_color_fea'][:2] for flake in flakes]) # Grayscale and V contrast
        
        # The average for each cluster gives the position for the center of that cluster (in the feature space)
        central_features = np.zeros([results['num_clusters'], consider_features.shape[1]])
        for i in range(results['num_clusters']):
             cluster_i = np.nonzero(assignment==i)[0]
             central_features[i,:] = np.mean(consider_features[cluster_i, :])
        to_save['sort_indices'] = np.argsort(np.abs(central_features).sum(1))
        to_save['unsort2sort'] = np.unique(to_save['sort_indices'], return_index=True)[1]
        
        to_save['cluster_centers'] = cluster_result.cluster_centers_[to_save['sort_indices']] # in (normed) feature space coordinates
        to_save['cluster_center_distances'] = cdist(to_save['cluster_centers'], to_save['cluster_centers'])
        
        
        
        # Output results
        ########################################
        
        # Save cluster results
        outfile = self.get_outfile(basename, output_dir, ext=run_args['file_extension'])
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'results of cluster analysis' ,
             'type' : 'data'
            } ,
            ]
        with open(outfile, 'wb') as fout:
            pickle.dump(to_save, fout)


        # Output images
        if run_args['verbosity']>=4:
            print('  Generating PCA 3D projection')

        # Pick a color for each cluster
        norm = mpl.colors.Normalize(vmin=0, vmax=results['num_clusters']-1)
        cmap = mpl.cm.jet
        #cmap = cmap_vge
        m = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        cluster_colors = [m.to_rgba(index) for index in range(results['num_clusters'])]
        
            
        pca = PCA(n_components=3)
        coordinates = pca.fit_transform(features)
        outfile = self.get_outfile(basename, output_dir, ext='-{}.png'.format(run_args['cluster_method']))
        self.plot_pca(outfile, coordinates, assignment, cluster_colors, **run_args)
        

        if run_args['verbosity']>=4:
            print('  Generating map of distances')
        
        outfile = self.get_outfile('distances', output_dir, ext='-{}.png'.format(run_args['cluster_method']))
        self.plot_distances(outfile, to_save['cluster_center_distances'], cluster_colors, **run_args)


        if run_args['verbosity']>=4:
            print('  Generating cluster images')
            
            
        self.plot_clusters(output_dir, to_save['cluster_centers'], **run_args)

        
        return results
    
    

    

    def plot_pca(self, outfile, coordinates, assignment, cluster_colors, **run_args):
        
        flake_colors = [cluster_colors[index] for index in assignment]

        # Centroid of each cluster (in PCA coordinates)
        num_clusters = np.max(assignment)+1
        cluster_coordinates = np.zeros([num_clusters, coordinates.shape[1]])
        for i in range(num_clusters):
             cluster_i = np.nonzero(assignment==i)[0]
             cluster_coordinates[i,:] = np.mean(coordinates[cluster_i,:], axis=0)
        cluster_index = range(cluster_coordinates.shape[0])
        
        
        plt.rcParams['xtick.labelsize'] = 8
        plt.rcParams['ytick.labelsize'] = 8
        plt.rcParams['axes.labelsize'] = 12
        plt.rcParams['lines.markersize'] = 5
        
        
        cmap = run_args['cmap'] if 'cmap' in run_args else 'jet'
        alpha = 0.12

        self.fig = plt.figure(figsize=(10,10))
        self.fig.subplots_adjust(left=0.08, right=0.95, bottom=0.08, top=0.95, hspace=0.15, wspace=0.15)
        
        
        self.ax = self.fig.add_subplot(2,2,2 , projection='3d')
        self.ax.scatter(coordinates[:,0], coordinates[:,1], coordinates[:,2], c=flake_colors, alpha=0.3)
        self.ax.set_xlabel('$\mathrm{PCA}_1$', labelpad=-4)
        self.ax.set_ylabel('$\mathrm{PCA}_2$', labelpad=-4)
        self.ax.set_zlabel('$\mathrm{PCA}_3$', labelpad=-2)
        self.ax.tick_params(axis='both', which='major', pad=-1)
        self.ax.view_init(elev=30, azim=45)
        
        
        self.ax = self.fig.add_subplot(2,2,1)
        self.ax.scatter(coordinates[:,0], coordinates[:,2], c=flake_colors, edgecolors=None, alpha=alpha)
        self.ax.set_xlabel('$\mathrm{PCA}_1$')
        self.ax.set_ylabel('$\mathrm{PCA}_3$')
        self.plot_cluster_number(cluster_coordinates, cluster_index, 0, 2, cluster_colors)
        xi, xf, yi, yf = self.ax.axis()
        self.ax.text(xi,yf, '{:,d} flakes in {} clusters'.format(len(assignment), len(cluster_colors)), size=10, verticalalignment='top', horizontalalignment='left', alpha=0.5)


        self.ax = self.fig.add_subplot(2,2,3)
        self.ax.scatter(coordinates[:,0], coordinates[:,1], c=flake_colors, cmap=cmap, edgecolors=None, alpha=alpha)
        self.ax.set_xlabel('$\mathrm{PCA}_1$')
        self.ax.set_ylabel('$\mathrm{PCA}_2$')
        self.plot_cluster_number(cluster_coordinates, cluster_index, 0, 1, cluster_colors)

        self.ax = self.fig.add_subplot(2,2,4)
        self.ax.scatter(coordinates[:,2], coordinates[:,1], c=flake_colors, cmap=cmap, edgecolors=None, alpha=alpha)
        self.ax.set_xlabel('$\mathrm{PCA}_3$')
        self.ax.set_ylabel('$\mathrm{PCA}_2$')
        self.plot_cluster_number(cluster_coordinates, cluster_index, 2, 1, cluster_colors)

        plt.savefig(outfile, dpi=300)
        plt.close()
        
    
    def plot_cluster_number(self, cluster_coordinates, cluster_index, coord1, coord2, cluster_colors):
        
        r = 0.3 # r=1 means no fade (strong color), r=0 means fully faded (appears white)
        cluster_colors_a = [ [ 1-(1-c[0])*r, 1-(1-c[1])*r, 1-(1-c[2])*r, c[3]] for c in cluster_colors]
        self.ax.scatter(cluster_coordinates[:,coord1], cluster_coordinates[:,coord2], s=25, c=cluster_colors_a, edgecolor=cluster_colors, alpha=1)
        
        for i in range(cluster_coordinates.shape[0]):
            self.ax.text(cluster_coordinates[i, coord1], cluster_coordinates[i, coord2], '{}'.format(i), size=3, horizontalalignment='center', verticalalignment='center')
        
        
    def plot_distances(self, outfile, cluster_center_distances, cluster_colors, plot_buffers=[0.15,0.05,0.15,0.05], **run_args):

        plt.rcParams['xtick.labelsize'] = 15
        plt.rcParams['ytick.labelsize'] = 15
        plt.rcParams['axes.labelsize'] = 20
        plt.rcParams['lines.markersize'] = 5


        self.fig = plt.figure( figsize=(8,8), facecolor='white' )
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )
        
        #plt.figtext(0,1, 'distances between clusters (in the feature space)', size=15, verticalalignment='top', horizontalalignment='left')
        self.ax.imshow(cluster_center_distances, cmap='viridis')
        
        self.ax.set_xlabel('cluster index')
        self.ax.set_ylabel('cluster index')


        xi, xf, yi, yf = self.ax.axis()
        
        s = 0.02
        n = len(cluster_colors)
        self.axt = self.fig.add_axes( [left_buf, bottom_buf+fig_height, fig_width, s] )
        self.axt.scatter(range(n), np.ones(n), c=cluster_colors)
        if n<160:
            for i in range(n):
                self.axt.text(i, 1, '{}'.format(i), size=4, horizontalalignment='center', verticalalignment='center')
        self.axt.axis([xi, xf, 0, 2])

        self.axt.axes.get_xaxis().set_visible(False)
        self.axt.axes.get_yaxis().set_visible(False)        
        
        self.axr = self.fig.add_axes( [left_buf+fig_width, bottom_buf, s, fig_height] )
        self.axr.scatter(np.ones(n), range(n), c=cluster_colors)
        if n<80:
            for i in range(n):
                self.axr.text(1, i, '{}'.format(i), size=4, horizontalalignment='center', verticalalignment='center')
        self.axr.axis([0, 2, yi, yf])
        
        self.axr.axes.get_xaxis().set_visible(False)
        self.axr.axes.get_yaxis().set_visible(False)
        
        plt.savefig(outfile, dpi=300)
    
        

    def plot_clusters(self, output_dir, cluster_centers, plot_buffers=[0.1,0.1,0.1,0.1], **run_args):

        plt.rcParams['xtick.labelsize'] = 15
        plt.rcParams['ytick.labelsize'] = 15
        plt.rcParams['axes.labelsize'] = 20
        plt.rcParams['lines.markersize'] = 5
        
        for i, center in enumerate(cluster_centers[:2]):
            
            if run_args['verbosity']>=5:
                print('    image for cluster {:03d}'.format(i))
            

            self.fig = plt.figure( figsize=(8,8), facecolor='white' )
            left_buf, right_buf, bottom_buf, top_buf = plot_buffers
            fig_width = 1.0-right_buf-left_buf
            fig_height = 1.0-top_buf-bottom_buf
            self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )            
            
            plt.figtext(0,1, 'cluster {:03d}'.format(i), size=20, verticalalignment='top', horizontalalignment='left')


            self.ax = self.fig.add_axes( [0, 0, 0.2, 0.95] )
            #self.ax.axes.get_xaxis().set_visible(False)
            #self.ax.axes.get_yaxis().set_visible(False)        
            self.ax.set_yticklabels([])
            self.ax.set_xticklabels([])
            vector = np.asarray([center]).transpose()
            self.ax.imshow(vector, cmap='inferno', aspect='auto')
            
            xi, xf, yi, yf = self.ax.axis()
            if len(center)<80:
                for ic, c in enumerate(center):
                    if c<0:
                        color = 'white'
                    else:
                        color = 'k'
                    self.ax.text((xi+xf)*0.5, ic, '{:.2f}'.format(c), color=color, size=8, verticalalignment='center', horizontalalignment='center')
            
            
            outfile = os.path.join(output_dir, 'cluster-{}-{:03d}.png'.format(run_args['cluster_method'], i))
            plt.savefig(outfile, dpi=300)
            plt.close(self.fig.number)

            
        
 
