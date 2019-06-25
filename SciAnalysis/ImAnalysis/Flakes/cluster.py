#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4

import pickle
from ..Protocols import *


from scipy.spatial.distance import cdist
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, MeanShift, estimate_bandwidth, AffinityPropagation, SpectralClustering # Clustering methods
from sklearn.decomposition import PCA

import skimage

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches

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
                        'feature_normed_range' : [-2, +2], # range for plotting the normed features
                        'bbox_pad' : 0.5,
                        'image_contrast' : (0, 1),
                        'image_contrast_trim' : None,
                        'overlays' : 3,
                        }
        self.run_args.update(kwargs)
        
        
        # WARNING: This association of features names is hard-coded, and is thus contingent
        # on the current implementation of Protocols.py>flake_analysis
        self.feature_names_color = [
            'g contrast',
            'v contrast',
            'gray',
            'gray std',
            'H',
            'S',
            'V',
            'H std',
            'S std',
            'V std',
            'R',
            'G',
            'B',
            'R std',
            'G std',
            'B std',
            'entropy'
            ]
        self.feature_names_color = self.feature_names_color + ['{}_inner'.format(f) for f in self.feature_names_color]
        
        self.feature_names_shape = ['P/A'] + ['hist {}'.format(i) for i in range(15)] + ['fractal dimension']
        
        
        
        
        
        
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
    
    def load_features(self, flakes, **run_args):
        
        if run_args['features']=='all':
            features = [ np.concatenate([flake['flake_color_fea'], flake['flake_shape_fea']]) for flake in flakes ]
            
            if 'flake_color_fea_names' in flakes[0]:
                self.feature_names_color = flakes[0]['flake_color_fea_names']
            if 'flake_shape_fea_names' in flakes[0]:
                self.feature_names_shape = flakes[0]['flake_shape_fea_names']
            
            self.feature_names = self.feature_names_color + self.feature_names_shape
            
        else:
            features = [ flake['flake_{}_fea'.format(run_args['features'])] for flake in flakes ]
            
            if run_args['features']=='color':
                if 'flake_color_fea_names' in flakes[0]:
                    self.feature_names = flakes[0]['flake_color_fea_names']
                else:
                    self.feature_names = self.feature_names_color
                    
            elif run_args['features']=='shape':
                if 'flake_shape_fea_names' in flakes[0]:
                    self.feature_names = flakes[0]['flake_shape_fea_names']
                else:
                    self.feature_names = self.feature_names_shape
                    
            else:
                if 'flake_{}_fea_names'.format(run_args['features']) in flakes[0]:
                    self.feature_names = flakes[0]['flake_{}_fea_names'.format(run_args['features'])]
                else:
                    self.feature_names = []
                    
        return np.asarray(features)


    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}
        clustering = {} # Save results of clustering operation
        
        # Aggregate results
        ########################################
        flakes = self.load_flakes(datas, **run_args)
        if run_args['verbosity']>=4:
            print('  {:,d} flakes identified in {:d} images'.format(len(flakes), len(datas)))
            
        features_orig = self.load_features(flakes, **run_args)
        
        
        
        # Clustering
        ########################################
        rescale = StandardScaler()
        features = rescale.fit_transform(features_orig)
        
        if run_args['verbosity']>=4:
            print("  Clustering {:,d} flakes using '{}'".format(len(flakes), run_args['cluster_method']))
            
        start = time.time()
        
        n_jobs = run_args['num_jobs'] if 'num_jobs' in run_args else -1
        
        if run_args['cluster_method']=='kmeans':
            cluster_result = KMeans(n_clusters=run_args['num_clusters'], random_state=0, n_jobs=n_jobs).fit(features)
            
        elif run_args['cluster_method']=='meanshift':
            bandwidth = estimate_bandwidth(features, quantile=0.1)#, n_samples=int(features.shape[0]/10))
            cluster_result = MeanShift(bandwidth=bandwidth, bin_seeding=True, n_jobs=n_jobs).fit(features)
            
        elif run_args['cluster_method']=='affinity':
            cluster_result = AffinityPropagation().fit(features)

        elif run_args['cluster_method']=='spectral':
            cluster_result = SpectralClustering(n_clusters=run_args['num_clusters'], n_jobs=n_jobs).fit(features)
            
        else:
            print("ERROR: clustering method '{}' not recognized.".format(run_args['cluster_method']))
            raise NotImplementedError

        clustering['cluster_result'] = cluster_result
        results['cluster_runtime'] = time.time()-start
        results['cluster_method'] = run_args['cluster_method']
            
        # Assignments are unsorted by default
        assignment = cluster_result.labels_
        results['num_clusters'] = len(np.unique(assignment))
        clustering['assignment'] = assignment # Label ids for each flake, saying what cluster it belongs to [unsorted indexing]

        if run_args['verbosity']>=4:
            print("    clustering took {:.1f}s ({:d} clusters)".format(results['cluster_runtime'], results['num_clusters']))
            
            
        # Sort clusters into a sensible order
        consider_features = np.asarray([flake['flake_color_fea'][:2] for flake in flakes]) # Grayscale and V contrast
        
        # The average for each cluster gives the position for the center of that cluster (in the feature space)
        central_features = np.zeros([results['num_clusters'], consider_features.shape[1]])
        for i in range(results['num_clusters']):
             cluster_i = np.nonzero(assignment==i)[0]
             central_features[i,:] = np.mean(consider_features[cluster_i, :])
        clustering['sort_indices'] = np.argsort(np.abs(central_features).sum(1))
        clustering['unsort2sort'] = np.unique(clustering['sort_indices'], return_index=True)[1]
        
        clustering['cluster_centers'] = cluster_result.cluster_centers_[clustering['sort_indices']] # in (normed) feature space coordinates [sorted indexing]
        clustering['cluster_centers_orig'] = rescale.inverse_transform(clustering['cluster_centers']) # in (original) feature space coordinates [sorted indexing]
        clustering['cluster_center_distances'] = cdist(clustering['cluster_centers'], clustering['cluster_centers']) # in (normed) feature space coordinates [sorted indexing]


        # Compute additional things
        ########################################
        # The distribution (histogram) for each feature dimension
        # Since these are normed they should look somewhat Gaussian
        vmin, vmax = run_args['feature_normed_range']
        distributions, dist_bin_edges = np.apply_along_axis(lambda x: np.histogram(x, bins=50, range=[vmin,vmax], density=True), 0, features)
        clustering['distributions'] = distributions
        clustering['dist_bin_edges'] = dist_bin_edges
        
        
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
            pickle.dump(clustering, fout)


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
        self.plot_distances(outfile, clustering['cluster_center_distances'], cluster_colors, **run_args)


        if run_args['verbosity']>=4:
            print('  Generating cluster images')
            
            
        self.plot_clusters(output_dir, clustering['cluster_centers'], clustering['cluster_centers_orig'], clustering['sort_indices'], distributions, dist_bin_edges, flakes, features, assignment, **run_args)

        
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
        self.ax.view_init(elev=30, azim=45-90)
        
        
        self.ax = self.fig.add_subplot(2,2,1)
        self.ax.scatter(coordinates[:,0], coordinates[:,2], c=flake_colors, edgecolors=None, alpha=alpha)
        self.ax.set_xlabel('$\mathrm{PCA}_1$')
        self.ax.set_ylabel('$\mathrm{PCA}_3$')
        self.overlay_cluster_number(cluster_coordinates, cluster_index, 0, 2, cluster_colors)
        xi, xf, yi, yf = self.ax.axis()
        self.ax.text(xi,yf, '{:,d} flakes in {} clusters'.format(len(assignment), len(cluster_colors)), size=10, verticalalignment='top', horizontalalignment='left', alpha=0.5)


        self.ax = self.fig.add_subplot(2,2,3)
        self.ax.scatter(coordinates[:,0], coordinates[:,1], c=flake_colors, cmap=cmap, edgecolors=None, alpha=alpha)
        self.ax.set_xlabel('$\mathrm{PCA}_1$')
        self.ax.set_ylabel('$\mathrm{PCA}_2$')
        self.overlay_cluster_number(cluster_coordinates, cluster_index, 0, 1, cluster_colors)

        self.ax = self.fig.add_subplot(2,2,4)
        self.ax.scatter(coordinates[:,2], coordinates[:,1], c=flake_colors, cmap=cmap, edgecolors=None, alpha=alpha)
        self.ax.set_xlabel('$\mathrm{PCA}_3$')
        self.ax.set_ylabel('$\mathrm{PCA}_2$')
        self.overlay_cluster_number(cluster_coordinates, cluster_index, 2, 1, cluster_colors)

        plt.savefig(outfile, dpi=300)
        plt.close()
        
    
    def overlay_cluster_number(self, cluster_coordinates, cluster_index, coord1, coord2, cluster_colors):
        
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
    
        

    def plot_clusters(self, output_dir, cluster_centers, cluster_centers_orig, sort_indices, distributions, dist_bin_edges, flakes, flake_features, assignment, plot_buffers=[0.01,0.0,0.0,0.045], **run_args):

        plt.rcParams['xtick.labelsize'] = 15
        plt.rcParams['ytick.labelsize'] = 15
        plt.rcParams['axes.labelsize'] = 20
        plt.rcParams['lines.markersize'] = 5
        
        #for i, feature_vector in enumerate(cluster_centers[:1]): # for testing
        for i, feature_vector in enumerate(cluster_centers):
            
            # i # [sorted indexing]
            # feature_vector # in (normed) feature space coordinates [sorted indexing]
            feature_vector_orig = cluster_centers_orig[i] # in (original) feature space coordinates [sorted indexing]
            
            i_before_sort = sort_indices[i] # [unsorted indexing]
            
            cluster_i = np.nonzero(assignment==i_before_sort)[0] # indices [in unsorted indexing] of all flakes matching this cluster
            flakes_cluster = np.asarray(flakes)[cluster_i] # flakes matching this cluster
            features_cluster = flake_features[cluster_i] # feature vectors matching this cluster
            
            
            self.plot_cluster(output_dir, '{:03d}'.format(i), feature_vector, feature_vector_orig, flakes_cluster, features_cluster, distributions, dist_bin_edges, plot_buffers=plot_buffers, **run_args)
            
            
    def plot_cluster(self, output_dir, cluster_name, feature_vector, feature_vector_orig, flakes_cluster, features_cluster, distributions, dist_bin_edges, plot_buffers=[0.01,0.0,0.0,0.045], **run_args):
                     
        if run_args['verbosity']>=5:
            print('    image for cluster {} ({:,d} flakes)'.format(cluster_name, len(flakes_cluster)))
        

        self.fig = plt.figure( figsize=(8,8), facecolor='white' )
        vmin, vmax = run_args['feature_normed_range']
        
        plt.figtext(0,1, 'cluster {} ({:,d} flakes)'.format(cluster_name, len(flakes_cluster)), size=20, verticalalignment='top', horizontalalignment='left')


        # Sidebar that shows the feature vector for the centroid of this cluster
        ########################################
        fea_w = 0.04
        fea_h = 0.95
        self.ax = self.fig.add_axes( [0.0, 0, fea_w, fea_h] )
        
        
        vector = np.asarray([feature_vector]).transpose()
        self.ax.imshow(vector, cmap='inferno', aspect='auto', vmin=vmin, vmax=vmax)
        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])
        
        xi, xf, yi, yf = self.ax.axis()
        if len(feature_vector)<80:
            for ifea, fea in enumerate(feature_vector):
                if fea<0:
                    color = 'white'
                else:
                    color = 'k'
                self.ax.text((xi+xf)*0.5, ifea, '{:.2f}'.format(fea), color=color, size=8, verticalalignment='center', horizontalalignment='center')
                
                self.ax.text(xf, ifea, '{:.3g}'.format(feature_vector_orig[ifea]), size=6, verticalalignment='center', horizontalalignment='left')
                
                # Miniature histogram (of the entire distribution)
                axc = self.fig.add_axes( [fea_w*2.2, fea_h-(ifea+1)*fea_h/len(feature_vector), fea_w*2.3, fea_h/len(feature_vector)] )
                w = dist_bin_edges[ifea][1]-dist_bin_edges[ifea][0]
                axc.bar( dist_bin_edges[ifea][:-1]+0.5*w, distributions[ifea], width=w, color='b', alpha=0.3 )
                plt.xlim(vmin,vmax)
                
                
                # Overlay the histogram for this cluster
                distribution, dist_bin_edge = np.histogram(features_cluster[:,ifea], bins=50, range=[vmin,vmax], density=True)
                distribution *= np.max(distributions[ifea])/np.max(distribution)
                #axc.bar( dist_bin_edge[:-1]+0.5*w, distribution, width=w, color='purple', alpha=0.2 )
                axc.plot( dist_bin_edge[:-1]+0.5*w, distribution, '-', color='purple', linewidth=0.8, alpha=0.3 )

                
                axc.axvline(fea, color='purple', linewidth=1)
                if fea<vmin:
                    axc.axvline(vmin, color='purple', linewidth=4)
                elif fea>vmax:
                    axc.axvline(vmax, color='purple', linewidth=4)
                
                axc.axes.get_xaxis().set_visible(False)
                axc.axes.get_yaxis().set_visible(False)

                if len(self.feature_names)==len(feature_vector):
                    axc.text(vmin, np.max(distributions[ifea]), self.feature_names[ifea], size=4, verticalalignment='top', horizontalalignment='left', alpha=0.25)
                axc.text(vmax, np.max(distributions[ifea]), '{:d}'.format(ifea), size=4, verticalalignment='top', horizontalalignment='right', alpha=0.25)
                
                
        
        # Images of example flakes for this cluster
        ########################################
        
        # Sort flakes by their distance from the cluster centroid (which is located at position "feature_vector")
        distances = cdist(features_cluster, [feature_vector], metric='euclidean')[:,0]
        sort_indices = np.argsort(distances)
        #time.sleep(100)
        flakes_cluster = flakes_cluster[sort_indices]
        #features_cluster = flakes_cluster[sort_indices]
        distances = distances[sort_indices]
        
        
        
        
        # The total area we have available for plotting flakes
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        left_buf += fea_w*( 2.2 + 2.3 )
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        #self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )
        
        
        
        
        nrows, ncols = 3, 7
        w = fig_width/ncols
        idx = 0
        ystart = bottom_buf+fig_height
        plt.figtext(left_buf, ystart, 'central', size=8, verticalalignment='bottom', horizontalalignment='left')
        for irow in range(nrows):
            for icol in range(ncols):
                ax_pos = [left_buf+icol*w, ystart-(irow+1)*w, w, w]
                self.plot_flake_image(ax_pos, flakes_cluster[idx], distances[idx], **run_args)
                idx += 1

        ystart = ystart-nrows*w - 0.015
        #nrows, ncols = 2, 6
        w = fig_width/ncols
        idx = int( np.clip( len(flakes_cluster)/2, idx, len(flakes_cluster)-nrows*ncols ) )
        plt.figtext(left_buf, ystart, 'generic', size=8, verticalalignment='bottom', horizontalalignment='left')
        for irow in range(nrows):
            for icol in range(ncols):
                ax_pos = [left_buf+icol*w, ystart-(irow+1)*w, w, w]
                self.plot_flake_image(ax_pos, flakes_cluster[idx], distances[idx], **run_args)
                idx += 1

        ystart = ystart-nrows*w - 0.015
        nrows, ncols = 2, 7
        w = fig_width/ncols
        idx = len(flakes_cluster)-nrows*ncols
        plt.figtext(left_buf, ystart, 'peripheral', size=8, verticalalignment='bottom', horizontalalignment='left')
        for irow in range(nrows):
            for icol in range(ncols):
                ax_pos = [left_buf+icol*w, ystart-(irow+1)*w, w, w]
                self.plot_flake_image(ax_pos, flakes_cluster[idx], distances[idx], **run_args)
                idx += 1
        
        
        outfile = os.path.join(output_dir, 'cluster-{}-{}.png'.format(run_args['cluster_method'], cluster_name))
        plt.savefig(outfile, dpi=300)
        plt.close(self.fig.number)

            
    def plot_flake_image(self, ax_pos, flake_i, distance, **run_args):
        
        # Load parent image
        img = plt.imread(flake_i['infile'])
        h, w, c = img.shape

        # Define image sub-region that has the flake in it
        y1, y2, x1, x2 = flake_i['bbox']
        # Make the crop border a bit bigger than the flake bounding box
        box_size = (1+run_args['bbox_pad'])*max( abs(x2-x1), abs(y2-y1) )
        x1p = int(np.clip((x1+x2)*0.5 - box_size/2, 0, w))
        x2p = int(np.clip((x1+x2)*0.5 + box_size/2, 0, w))
        y1p = int(np.clip((y1+y2)*0.5 - box_size/2, 0, h))
        y2p = int(np.clip((y1+y2)*0.5 + box_size/2, 0, h))
        box = y1p, y2p, x1p, x2p

        # Adjust image of flake
        flake = img[y1p:y2p , x1p:x2p, :]
        in_range = self.get_in_range(img, run_args['image_contrast'], **run_args)
        flake = skimage.exposure.rescale_intensity(flake, in_range=in_range, out_range='dtype')
        
        
        # Plot flake
        self.ax = self.fig.add_axes(ax_pos)
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)
        
        self.ax.imshow(flake)

        xi, xf, yi, yf = self.ax.axis()
        yc, xc = flake_i['center_of_mass']
        s = '{}\nflake{:03d}\n({}, {})'.format(flake_i['infile'], flake_i['index'], int(xc), int(yc))
        self.ax.text(xi, yf, s, color='white', size=3, verticalalignment='top', horizontalalignment='left')
        

        self.ax.text(xi, yi, '${:.1f} \, \mathrm{{\mu m}}$'.format(flake_i['radius_um']), color='r', size=5, verticalalignment='bottom', horizontalalignment='left')
        
        self.ax.text((xi+xf)*0.5, yi, '{:.1f}'.format(distance), color='white', size=2, verticalalignment='bottom', horizontalalignment='center')

        self.ax.text(xf, yi, '{:.3f}'.format(flake_i['flake_contrast']), color='orange', size=3, verticalalignment='bottom', horizontalalignment='right')

        # Various overlays on the flake
        xc -= x1p
        yc -= y1p
        size = flake_i['radius_pixels']
        
        if run_args['overlays']>=1:
            c = flake_i['contour']
            xs = (c[:,0] - x1p)
            ys = (c[:,1] - y1p)
            self.ax.plot(xs, ys, '-', linewidth=0.6, color='r', dashes=[4,1], alpha=0.2)
        
        if run_args['overlays']>=7:
            c = flake_i['convex_hull']
            xs = (c[:,1] - x1p)
            ys = (c[:,0] - y1p)
            self.ax.plot(xs, ys, '-', linewidth=0.5, color='g', alpha=0.5)

        if run_args['overlays']>=5:
            rect = patches.Rectangle( ((x1-x1p), (y1-y1p)), (x2-x1), (y2-y1), linewidth=1.0, edgecolor='orange', facecolor='none', alpha=0.5)
            self.ax.add_patch(rect)
        

        if run_args['overlays']>=3:
            # Cross hair and circle
            rect = patches.Rectangle( (xc-size/2, yc), size, 0, linewidth=0.6, edgecolor='r', facecolor='none', alpha=0.3) # Horizontal bar
            self.ax.add_patch(rect)
            rect = patches.Rectangle( (xc, yc-size/2), 0, size, linewidth=0.6, edgecolor='r', facecolor='none', alpha=0.3) # Vertical bar
            self.ax.add_patch(rect)

        if run_args['overlays']>=5:
            # Circle overlay denoting size
            circ = patches.Circle(xy=(xc,yc), radius=size, linewidth=0.6, edgecolor='r', facecolor='none', alpha=0.3)
            self.ax.add_patch(circ)
        

        
    def get_in_range(self, data, im_contrast, image_contrast_trim=None, **run_args):
        
        if image_contrast_trim is not None:
            image_contrast_trim = np.clip(image_contrast_trim, 0, 0.95)
            
            avg = np.average(data)
            avg /= 255
            amt = image_contrast_trim
            im_contrast = ( avg*amt , 1.0-(1.0-avg)*amt )

        in_range = ( im_contrast[0]*255, im_contrast[1]*255 )    
        
        return in_range




class select_flakes(cluster):
    
    def __init__(self, name='select_flakes', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.pkl'
        self.run_args = {
                        'file_extension' : '.pkl',
                        'force' : False,
                        'verbosity' : 3,
                        'num_jobs' : None,
                        'num_clusters' : 20,
                        'cluster_method' : 'selection',
                        'features' : 'all', # 'shape', 'color', 'all'
                        'feature_normed_range' : [-2, +2], # range for plotting the normed features
                        'bbox_pad' : 0.5,
                        'image_contrast' : (0, 1),
                        'image_contrast_trim' : None,
                        'overlays' : 3,
                        }
        self.run_args.update(kwargs)
    
    
        # WARNING: This association of features names is hard-coded, and is thus contingent
        # on the current implementation of Protocols.py>flake_analysis
        self.feature_names_color = [
            'g contrast',
            'v contrast',
            'gray',
            'gray std',
            'H',
            'S',
            'V',
            'H std',
            'S std',
            'V std',
            'R',
            'G',
            'B',
            'R std',
            'G std',
            'B std',
            'entropy'
            ]
        self.feature_names_color = self.feature_names_color + ['{}_inner'.format(f) for f in self.feature_names_color]
        
        self.feature_names_shape = ['P/A'] + ['hist {}'.format(i) for i in range(15)] + ['fractal dimension']

    
    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}

        # Load all flakes identified by "find_flakes" protocol
        flakes = self.load_flakes(datas, **run_args)
        if run_args['verbosity']>=4:
            print('  {:,d} flakes identified in {:d} images'.format(len(flakes), len(datas)))


        # Load data aggregated from the "cluster" protocol
        savefile = self.get_outfile(basename, output_dir+'/../cluster/', ext=run_args['file_extension'])
        with open(savefile, 'rb') as fin:
            clustering = pickle.load(fin)


        # Compute the rescaling of feature vectors
        features_all_orig = self.load_features(flakes, **run_args)
        rescale = StandardScaler()
        features_all_rescaled = rescale.fit_transform(features_all_orig)

        flakes_selected = self.select(flakes, features_all_orig, features_all_rescaled, **run_args)
        features_orig = self.load_features(flakes_selected, **run_args)
        
        
        feature_vector_orig = np.average(features_orig, axis=0)
        feature_vector = rescale.transform( [feature_vector_orig] )[0]
        features = rescale.transform(features_orig)
        
        
        if run_args['verbosity']>=4:
            print("  Selected {:,d} flakes using '{}'".format(len(flakes_selected), run_args['cluster_method']))
        
        
        self.plot_cluster(output_dir, cluster_name='target', feature_vector=feature_vector, feature_vector_orig=feature_vector_orig, flakes_cluster=flakes_selected, features_cluster=features, distributions=clustering['distributions'], dist_bin_edges=clustering['dist_bin_edges'], **run_args)


        return results


    def extract_feature(self, feature_name, flake_features, **run_args):

        i = self.feature_names.index(feature_name)
        if run_args['verbosity']>=5:
            print("    Extracting {}, index {}".format(feature_name, i))
                  
        return flake_features[:,i]
                    

    def select(self, flakes, flake_features_orig, flake_features_rescaled, **run_args):
        
        radius_um = np.asarray( [ f['radius_um'] for f in flakes ] )
        flake_contrast = np.asarray( [ f['flake_contrast'] for f in flakes ] )
        
        g_std = self.extract_feature('gray std', flake_features_rescaled, **run_args)
        g = self.extract_feature('gray', flake_features_rescaled, **run_args)
        g_std_rel = g_std/g
        
        
        
        conditions = [
                    (radius_um>2.0) ,
                    (radius_um<40.0) ,
                    (flake_contrast>-0.02) ,
                    (flake_contrast<-0.01) ,
                    (g_std_rel<0.1) ,
                    
                    ]
            
        
        #idx = np.where( (radius_um>5.0) & (radius_um<100.0) & (flake_contrast>-0.2) & (flake_contrast<-0.01) )[0]
        idx = np.where(np.all(conditions, axis=0))[0]
        
        flakes = np.asarray(flakes)[idx]
        
        return flakes



            
