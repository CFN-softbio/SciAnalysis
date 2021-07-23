#!/usr/bin/python3
import numpy as np


class CoordinateTransformation():
    
    def __init__(self, para1, para2, x_range=[-2, 52], y_range=[-2, 52], name='CoordTrans', method='fsolve'):
        
        self.para1, self.para2 = para1, para2
        self.x_range, self.y_range = x_range, y_range
        self.name = name
        self.prefix = '{}> '.format(self.name)
        self.method = method
        
        self.lookup_matrices = None
        self.linear_h2x, self.linear_h1y = None, None
        self.linear_h1_model, self.linear_h2_model = None, None
        if self.method=='lookup':
            self.generate_translation_maps()
        elif self.method=='linear':
            self.generate_linear_fit()
        elif self.method=='linear_model':
            self.generate_linear_model()
            
        self.determine_bounds()
            

    def val_stats(self, values, name='z'):
        '''Helper function to output the value ranges of an array.'''
        span = np.max(values)-np.min(values)
        print("  {} = {:.2g} ± {:.2g} (span {:.2g}, from {:.3g} to {:.3g})".format(name, np.average(values), np.std(values), span, np.min(values), np.max(values)))                
        
            
    def generate_translation_maps(self, h_step=0.1, d_step=0.1, d_extra=28):
        '''Generate maps in thickness space.'''

        xi, xf = self.x_range[0]-d_extra, self.x_range[1]+d_extra
        yi, yf = self.y_range[0]-d_extra, self.y_range[1]+d_extra
        xv, yv = np.arange(xi, xf, d_step), np.arange(yi, yf, d_step)
        X, Y = np.meshgrid(xv, yv)
        hL1, hL2 = self.xy_to_hh(X, Y)
        
        #self.val_stats(hL1, 'h1')
        #self.val_stats(hL2, 'h2')

        # The thickness space we will project into
        h1v = np.arange(np.min(hL1)+h_step, np.max(hL1)-h_step, h_step)
        h2v = np.arange(np.min(hL2)+h_step, np.max(hL2)-h_step, h_step)
        H1, H2 = np.meshgrid(h1v, h2v)
        #val_stats(H1, 'H1')
        #val_stats(H2, 'H2')
        
        
        print('{}Generating interpolation maps ({:d}×{:d} = {:,d} points)...'.format(self.prefix, H1.shape[0], H1.shape[1], H1.size))
        # Generate maps, which we will use as lookup values
        # Conceptually, this is like a surface plot, where the two axes are
        # h1 and h2, and the height of the surface is the corresponding x value
        # (or y value)
        import scipy.interpolate
        import time
        start_time = time.time()
        
        POINTS = np.column_stack((hL1.ravel(), hL2.ravel()))
        VALUES = X.ravel()
        X_MAP = scipy.interpolate.griddata(POINTS, VALUES, (H1, H2), method='linear')
        #X_MAP = np.ma.masked_where( np.isnan(X_MAP), X_MAP)
        #val_stats(X_MAP, 'X_MAP')

        VALUES = Y.ravel()
        Y_MAP = scipy.interpolate.griddata(POINTS, VALUES, (H1, H2), method='linear')
        #Y_MAP = np.ma.masked_where( np.isnan(Y_MAP), Y_MAP)
        #val_stats(Y_MAP, 'Y_MAP')
        
        
        took = time.time()-start_time
        print('{}    took {:.1f} s ({:.2f} s/pt)'.format(self.prefix, took, took/Y_MAP.size))
        
        self.lookup_matrices = h1v, h2v, X_MAP, Y_MAP


    def get_lookup_matrices(self):
        if self.lookup_matrices is None:
            self.generate_translation_maps()
            
        return self.lookup_matrices
        
        
    def xy_to_hh(self, x_vals, y_vals):
        '''Convert from sample (x,y) coordinates into calibrated
        film thickness (h1,h2) values.'''
        
        a1, b1, c1, d1 = self.para1
        a2, b2, c2, d2 = self.para2
        
        h_L1 = a1*np.square(y_vals) - b1*y_vals + c1*x_vals + d1
        h_L2 = a2*np.square(x_vals) - b2*x_vals + c2*y_vals + d2

        #h_total = h_L1 + h_L2
        #phi_C67 = h_L2/h_total

        return h_L1, h_L2
    
    
    def hh_to_xy(self, h1, h2):
        '''Convert from thickness (h1,h2) into sample position (x,y).
        We do this using a lookup table.'''
        
        compute_method = getattr(self, f'hh_to_xy_{self.method}')
        #if self.method=='lookup':
            #x_predicted, y_predicted = self.hh_to_xy_lookup(h1, h2)
        x_predicted, y_predicted = compute_method(h1, h2)
            
        return x_predicted, y_predicted

    def hh_to_xy_array(self, h1, h2):
        xs = np.zeros_like(h1)
        ys = np.zeros_like(h1)
        for i, (h1c, h2c) in enumerate(zip(h1, h2)):
            xs[i], ys[i] = self.hh_to_xy(h1c, h2c)
            
        return xs, ys

    def hh_to_xy_2Darray(self, H1, H2, h1_axis, h2_axis):
        X = np.zeros_like(H2)
        Y = np.zeros_like(H1)
        for i, h1 in enumerate(h1_axis):
            h1_cur = np.ones(len(h2_axis))*h1
            x, y = self.hh_to_xy_array(h1_cur, h2_axis)
            X[i] = x
            Y[i] = y
            
        return X, Y

    
    
    def generate_linear_fit(self, num=50):
        from scipy.stats import linregress

        # Estimate how h_L2 varies with x
        x = np.linspace(self.x_range[0], self.x_range[1], num=num, endpoint=True)
        y = np.ones_like(x)*(self.y_range[0]+self.y_range[1])*0.5 # Assume center of range
        h_L1, h_L2 = self.xy_to_hh(x, y)
        slope, intercept, r_value, p_value, std_err = linregress(x,h_L2)
        self.linear_h2x = [slope, intercept]

        # Estimate how h_L1 varies with y
        y = np.linspace(self.y_range[0], self.y_range[1], num=num, endpoint=True)
        x = np.ones_like(y)*(self.x_range[0]+self.x_range[1])*0.5 # Assume center of range
        h_L1, h_L2 = self.xy_to_hh(x, y)
        slope, intercept, r_value, p_value, std_err = linregress(y,h_L1)
        
        self.linear_h1y = [slope, intercept]


    def get_linear_fit(self, num=50):
        if self.linear_h2x is None or self.linear_h1y is None:
            self.generate_linear_fit(num=num)
            
        return self.linear_h2x, self.linear_h1y

    def generate_linear_model(self, num=50):
        # Build a multi-dimensional linear model within the region of interest
        
        x = np.linspace(self.x_range[0], self.x_range[1], num=num, endpoint=True)
        y = np.linspace(self.y_range[0], self.y_range[1], num=num, endpoint=True)
        X, Y = np.meshgrid(x, y)
        
        H1, H2 = self.xy_to_hh(X, Y)
        
        from sklearn import linear_model
        
        XY = np.column_stack( (X.ravel(), Y.ravel()) ) # Shape suitable for regression

        reg = linear_model.LinearRegression()
        reg.fit( XY, H1.ravel() )
        #ms, b = reg.coef_, reg.intercept_
        #predict_h1 = reg.predict( np.asarray([x,y]).reshape(1,-1) )
        #predict_h1 = ms[0]*x + ms[1]*y + b
        self.linear_h1_model = reg

        reg = linear_model.LinearRegression()
        reg.fit( XY, H2.ravel() )
        self.linear_h2_model = reg

    def get_linear_model(self, num=50):
        if self.linear_h1_model is None or self.linear_h2_model is None:
            self.generate_linear_model(num=num)
            
        return self.linear_h1_model, self.linear_h2_model

    
    def hh_to_xy_linear(self, h1, h2, num=50):
        
        linear_h2x, linear_h1y = self.get_linear_fit(num=num)
        
        # y = m*x + b --> h_L2 = slope*x + intercept
        slope, intercept = linear_h2x
        x_predicted = (h2-intercept)/slope
        
        # y = m*x + b --> h_L1 = slope*y + intercept
        slope, intercept = linear_h1y
        y_predicted = (h2-intercept)/slope
        
        return x_predicted, y_predicted


    def hh_to_xy_linear_model(self, h1, h2, num=50):
        
        linear_h1_model, linear_h2_model = self.get_linear_model(num=num)
        (m_x1, m_y1), b_1 = linear_h1_model.coef_, linear_h1_model.intercept_
        (m_x2, m_y2), b_2 = linear_h2_model.coef_, linear_h2_model.intercept_
        
        # h1 = m_x1*x + m_y1*y + b_1
        # h2 = m_x2*x + m_y2*y + b_2
        
        # Find the intersection of the two "lines" in (x,y) space that
        # satsify the linear_model. This intersection point is the
        # solution that satisfies both h1 and h2 models.
        x = ( m_y1*(h2-b_2) - m_y2*(h1-b_1) )/( m_y1*m_x2 - m_y2*m_x1 )
        y = (h1 - m_x1*x - b_1)/m_y1
        #y = (h2 - m_x2*x - b_2)/m_y2 # Either equation will work
        
        return x, y
    
    
    def hh_to_xy_fsolve(self, h1, h2, guess='linear'):
        
        if guess=='center':
            x_guess = (self.x_range[0]+self.x_range[1])*0.5
            y_guess = (self.y_range[0]+self.y_range[1])*0.5
        elif guess=='low':
            x_guess, y_guess = self.x_range[1], self.y_range[0]
        elif guess=='high':
            x_guess, y_guess = self.x_range[0], self.y_range[1]
        else: # Default is 'linear'
            x_guess, y_guess = self.hh_to_xy_linear(h1, h2)
        
        from scipy.optimize import fsolve
        
        a1, b1, c1, d1 = self.para1
        a2, b2, c2, d2 = self.para2
        
        
        def func(x):
            v1 = a1*np.square(x[1]) - b1*x[1] + c1*x[0] + d1 - h1
            v2 = a2*np.square(x[0]) - b2*x[0] + c2*x[1] + d2 - h2
            return [v1, v2]
        
        x_predicted, y_predicted = fsolve(func, np.asarray([x_guess,y_guess]))
        
        return x_predicted, y_predicted


    def hh_to_xy_root(self, h1, h2):
        
        x_guess, y_guess = self.hh_to_xy_linear(h1, h2)
        
        from scipy.optimize import root
        
        a1, b1, c1, d1 = self.para1
        a2, b2, c2, d2 = self.para2
        
        
        def func(x):
            v1 = a1*np.square(x[1]) - b1*x[1] + c1*x[0] + d1 - h1
            v2 = a2*np.square(x[0]) - b2*x[0] + c2*x[1] + d2 - h2
            return [v1, v2]
        
        sol = root(func, np.asarray([x_guess,y_guess]))
        x_predicted, y_predicted = sol.x
        
        return x_predicted, y_predicted
    
        
    def hh_to_xy_lookup(self, h1, h2):
        
        lookup_matrices = self.get_lookup_matrices()
        
        try:
            h1v, h2v, X_MAP, Y_MAP = lookup_matrices
            
            # Find indices in the map
            idx_h1 = np.where(h1v>=h1)[0][0]
            idx_h2 = np.where(h2v>=h2)[0][0]
            x_predicted, y_predicted = X_MAP[idx_h2,idx_h1], Y_MAP[idx_h2,idx_h1]
            
            if np.isnan(x_predicted) or np.isnan(y_predicted):
                print('      hh_to_xy produced nan')
                x_predicted, y_predicted = None, None
            
        except IndexError as e:
            print('      hh_to_xy raised exception: {}'.format(e.__class__.__name__))
            x_predicted, y_predicted = None, None
        
        
        return x_predicted, y_predicted


    def determine_bounds(self, num=20, method='detaulf'):
        # Figure out what the bounds should be in (h1,h2) space
        
        print('x bounds: {:.1f} to {:.1f}'.format(self.x_range[0], self.x_range[1]))
        print('y bounds: {:.1f} to {:.1f}'.format(self.y_range[0], self.y_range[1]))
        
        
        if method=='corner':
            x = np.linspace(self.x_range[0], self.x_range[1], num=num, endpoint=True)
            y = np.linspace(self.y_range[0], self.y_range[1], num=num, endpoint=True)
            h1, h2 = self.xy_to_hh(x, y) # Diagonal across bounds
            h1d, h2d = self.xy_to_hh(x, y[::-1]) # Other diagonal
            
            h1_min, h1_max = np.min([h1, h1d]), np.max([h1, h1d])
            h2_min, h2_max = np.min([h2, h2d]), np.max([h2, h2d])
            
            self.h1_range, self.h2_range = [h1_min, h1_max], [h2_min, h2_max]
            
        else:
            # default: generate a grid and search for extrema
            h1_min, h1_max = 1e9, 0
            h2_min, h2_max = 1e9, 0
            for x in np.linspace(self.x_range[0], self.x_range[1], num=num, endpoint=True):
                for y in np.linspace(self.y_range[0], self.y_range[1], num=num, endpoint=True):
                    h1, h2 = self.xy_to_hh(x, y)
                    x_check, y_check = self.hh_to_xy(h1, h2)
                    #if x_check is None or y_check is None:
                        #print('xy({:.1f}, {:.1f}) -> hh[{:.1f},{:.1f}] -> xy({},{})'.format(x, y, h1, h2, x_check, y_check))
                    #else:
                        #print('xy({:.1f}, {:.1f}) -> hh[{:.1f},{:.1f}] -> xy({:.1f},{:.1f})'.format(x, y, h1, h2, x_check, y_check))
                    h1_min = min(h1_min, h1)
                    h1_max = max(h1_max, h1)
                    h2_min = min(h2_min, h2)
                    h2_max = max(h2_max, h2)
            
            self.h1_range, self.h2_range = [h1_min, h1_max], [h2_min, h2_max]
            
        print('h1 bounds: {:.1f} to {:.1f}'.format(self.h1_range[0], self.h1_range[1]))
        print('h2 bounds: {:.1f} to {:.1f}'.format(self.h2_range[0], self.h2_range[1]))



        
    def in_xy_bounds(self, x, y):
        if x<self.x_range[0] or y<self.y_range[0] or x>self.x_range[1] or y>self.y_range[1]:
            return False
        else:
            return True

    def in_hh_box(self, h1, h2):
        if h1<self.h1_range[0] or h2<self.h2_range[0] or h1>self.h1_range[1] or h2>self.h2_range[1]:
            return False
        else:
            return True
            
    def in_hh_polygon(self, h1, h2):
        x, y = self.hh_to_xy(h1, h2)
        return self.in_xy_bounds(x, y)
    
    def in_hh_valid(self, h1, h2, threshold=1, **kwargs):
        if not self.in_hh_box(h1, h2):
            return False
        if not self.in_hh_polygon(h1, h2):
            return False
        distortion = self.distortion(h1, h2, **kwargs)
        return distortion<threshold
        
        

    def confirm_working(self, h_range=None, num=4):
        self.determine_bounds()
        
        if h_range is None:
            h1_range, h2_range = self.h1_range, self.h2_range
        else:
            h1_range, h2_range = h_range, h_range
        
        deltas = []
        
        print('{}Checking transformations for {:d}×{:d} = {:,d} points...'.format(self.prefix, num, num, num*num))
        import time
        start_time = time.time()
        
        # Confirm that maps works:
        for h1_try in np.linspace(h1_range[0], h1_range[1], num=num):
            for h2_try in np.linspace(h2_range[0], h2_range[1], num=num):
                
                x_predicted, y_predicted = self.hh_to_xy(h1_try, h2_try)
                
                if x_predicted is not None and y_predicted is not None:
                    h1_check, h2_check = self.xy_to_hh(x_predicted, y_predicted)
                    
                    del1, del2 = h1_check-h1_try, h1_check-h1_try


                    b = '' if self.in_xy_bounds(x_predicted, y_predicted) else '    out-of-bounds'
                    print('    Point hh[{:.1f},{:.1f}] -> xy({:.1f},{:.1f}) -> hh[{:.1f},{:.1f}]    [{:+1.2f},{:+1.2f}]{}'.format(h1_try, h2_try, x_predicted, y_predicted, h1_check, h2_check, del1, del2, b))
                    
                    deltas.append(abs(del1))
                    deltas.append(abs(del2))
                        
                else:
                    print('    Point hh[{:.1f},{:.1f}] does not work'.format(h1_try, h2_try))
                    
        self.val_stats(np.asarray(deltas), 'differences')

        took = time.time()-start_time
        print('{}    took {:.3f} s ({:.3f} s/pt)'.format(self.prefix, took, took/(num*num)))


    def plot_maps(self, name='sample', gnum=4, num=1000, d_extra=20):
        # Plot the individual conversion equations
        import matplotlib as mpl
        mpl.rcParams['mathtext.fontset'] = 'cm'
        import pylab as plt
        
        self.fig = plt.figure(figsize=(10,5))
        plot_buffers = [0.1,0.01,0.1,0.05]
        gap = 0.02
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 0.5-right_buf-left_buf-gap
        fig_height = 1.0-top_buf-bottom_buf
        self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )
        self.ax2 = self.fig.add_axes( [2*left_buf+fig_width+gap, bottom_buf, fig_width, fig_height] )
        

        conditions = [
            ['k', 2.0, gnum] , # Vertical major lines
            ['k', 0.5, gnum*4] , # Vertical minor lines
            ]
        for c, t, cnum in conditions:
            for x in np.linspace(self.x_range[0], self.x_range[1], num=cnum, endpoint=True):
                y_plot = np.linspace(self.y_range[0], (self.y_range[1]+self.y_range[0])/2, num=num, endpoint=True)
                x_plot = np.ones_like(y_plot)*x
                self.ax.plot(x_plot, y_plot, '-', color=c, linewidth=t, dashes=[2,1])
                h1, h2 = self.xy_to_hh(x_plot, y_plot)
                self.ax2.plot(h2, h1, '-', color=c, linewidth=t, dashes=[2,1])

                y_plot = np.linspace((self.y_range[1]+self.y_range[0])/2, self.y_range[1], num=num, endpoint=True)
                x_plot = np.ones_like(y_plot)*x
                self.ax.plot(x_plot, y_plot, '-', color=c, linewidth=t)
                h1, h2 = self.xy_to_hh(x_plot, y_plot)
                self.ax2.plot(h2, h1, '-', color=c, linewidth=t)


        conditions = [
            ['0.5', 2.0, gnum] , # Horizontal major lines
            ['0.5', 0.5, gnum*4] , # Horizontal minor lines
            ]
        for c, t, cnum in conditions:
            for y in np.linspace(self.y_range[0], self.y_range[1], num=cnum, endpoint=True):
                x_plot = np.linspace(self.x_range[0], (self.x_range[1]+self.x_range[0])/2, num=num, endpoint=True)
                y_plot = np.ones_like(x_plot)*y
                self.ax.plot(x_plot, y_plot, '-', color=c, linewidth=t)
                h1, h2 = self.xy_to_hh(x_plot, y_plot)
                self.ax2.plot(h2, h1, '-', color=c, linewidth=t)

                x_plot = np.linspace((self.x_range[1]+self.x_range[0])/2, self.x_range[1], num=num, endpoint=True)
                y_plot = np.ones_like(x_plot)*y
                self.ax.plot(x_plot, y_plot, '-', color=c, linewidth=t, dashes=[2,1])
                h1, h2 = self.xy_to_hh(x_plot, y_plot)
                self.ax2.plot(h2, h1, '-', color=c, linewidth=t, dashes=[2,1])
            
        
        conditions = [
            ['b', 2.0, gnum] ,
            ['b', 0.5, gnum*4] ,
            ]
        for c, t, cnum in conditions:
            for h1 in np.linspace(self.h1_range[0], self.h1_range[1], num=cnum, endpoint=True):
                h2_plot = np.linspace(self.h2_range[0], self.h2_range[1], num=num, endpoint=True)
                h1_plot = np.ones_like(h2_plot)*h1
                self.ax2.plot(h2_plot, h1_plot, '-', color=c, linewidth=t, alpha=0.2, zorder=-10)
                
                x, y = self.hh_to_xy_array(h1_plot, h2_plot)
                self.ax.plot(x, y, '-', color=c, linewidth=t, alpha=0.2, zorder=-10)

        conditions = [
            ['purple', 2.0, gnum] ,
            ['purple', 0.5, gnum*4] ,
            ]
        for c, t, cnum in conditions:
            for h2 in np.linspace(self.h2_range[0], self.h2_range[1], num=cnum, endpoint=True):
                h1_plot = np.linspace(self.h1_range[0], self.h1_range[1], num=num, endpoint=True)
                h2_plot = np.ones_like(h1_plot)*h2
                self.ax2.plot(h2_plot, h1_plot, '-', color=c, linewidth=t, alpha=0.2, zorder=-10)
                
                x, y = self.hh_to_xy_array(h1_plot, h2_plot)
                self.ax.plot(x, y, '-', color=c, linewidth=t, alpha=0.2, zorder=-10)


        self.ax.set_xlabel('$x \, (\mathrm{mm}$)', size=14)
        self.ax.set_ylabel('$y \, (\mathrm{mm}$)', size=14)


        self.ax2.set_xlabel('$h_{\mathrm{L2}} \, (\mathrm{nm})$', size=14)
        self.ax2.set_ylabel('$h_{\mathrm{L1}} \, (\mathrm{nm})$', size=14)        
        
        self.ax.set_aspect('equal')
        self.ax2.set_aspect('equal')
        
        self.fig.savefig('{}-maps.png'.format(name), dpi=300)
        plt.close()


    def distortion(self, h1, h2, n_power=1):
        '''Return the distortion of the data in this part of the coordinate transformation.
        The returned value is 0 if the data is undistorted (viable).
        A value >1 indicates unacceptable amounts of distortion.'''


        # Compute gradients
        from scipy.optimize import approx_fprime
        
        def hh_to_x(xk):
            return self.hh_to_xy(xk[0], xk[1])[0]
        def hh_to_y(xk):
            return self.hh_to_xy(xk[0], xk[1])[1]
        
        
        span = np.average( [ self.h1_range[1]-self.h1_range[0], self.h2_range[1]-self.h2_range[0] ] )
        epsilon = [span/50, span/50]
        
        gradX = approx_fprime((h1,h2), hh_to_x, epsilon)
        gradX = np.sqrt( np.square(gradX[0]) + np.square(gradX[1]) )
        gradY = approx_fprime((h1,h2), hh_to_y, epsilon)
        gradY = np.sqrt( np.square(gradY[0]) + np.square(gradY[1]) )
        gradT = np.sqrt( np.square(gradX) + np.square(gradY) )

        # Renormalize gradient to yield a signal for "how nonlinear is the transform?"
        (slope_h2x, intercept_h2x), (slope_h1y, intercept_h1y) = self.get_linear_fit()
        slope = np.average((np.abs(slope_h2x), np.abs(slope_h1y))) # nm/mm
        grad_typical = 1/slope # mm/nm
        gradN = gradT/grad_typical
        # gradN = 1 is normal/expected; anything
        # gradN > 1 represents high distortion/compression; h(x) was flat, x(h) is steep
        # gradN < 1 represents high distortion/stretching; h(x) was vertical, x(h) is flat
        gradN = (gradN if gradN>=1 else 1/gradN) - 1
        # Now gradN is [0..+inf]; gradN~1 means gradient was twice as large as expectations
        #gradN *= 0.5

        # Cyclic delta
        x, y = self.hh_to_xy(h1, h2)
        h1c, h2c = self.xy_to_hh(x, y)
        delta = np.sqrt( np.square(h1-h1c) + np.square(h2-h2c) )
        delta *= 100/span # I.e. we consider a deviation of 1/100 of the range to be 'large'

        # Difference from linear
        XY = np.column_stack( (x, y) )
        linear_h1_model, linear_h2_model = self.get_linear_model()
        #x, y = self.hh_to_xy_linear_model(h1, h2)

        predict_H1 = linear_h1_model.predict(XY)[0]
        predict_H2 = linear_h2_model.predict(XY)[0]
        diff_H1 = (h1-predict_H1)/span
        diff_H2 = (h2-predict_H2)/span

        diff = np.sqrt( np.square(diff_H1) + np.square(diff_H2) )*2

        d = np.power(gradN + delta + diff, n_power)
        
        return d


    def distortion_array(self, h1, h2, **kwargs):
        d = np.zeros_like(h1)
        for i, (h1c, h2c) in enumerate(zip(h1, h2)):
            d[i] = self.distortion(h1c, h2c, **kwargs)
            
        return d

    def distortion_2Darray(self, H1, H2, h1_axis, h2_axis, **kwargs):
        D = np.zeros_like(H1)
        for i, h1 in enumerate(h1_axis):
            h1_cur = np.ones(len(h2_axis))*h1
            d = self.distortion_array(h1_cur, h2_axis, **kwargs)
            D[i] = d
            
        return D

    def distortion_2Darray_map(self, H1, H2, h1_axis, h2_axis, n_power=1, **kwargs):
        
        X, Y = self.hh_to_xy_2Darray(H1, H2, h1_axis, h2_axis)
        
        dh1 = (h1_axis[-1]-h1_axis[0])/(len(h1_axis)-1)
        dh2 = (h2_axis[-1]-h2_axis[0])/(len(h2_axis)-1)
        
        # Gradients
        gradX = np.gradient(X)
        gradX = np.sqrt( np.square(gradX[0]/dh1) + np.square(gradX[1]/dh2) )
        gradY = np.gradient(Y)
        gradY = np.sqrt( np.square(gradY[0]/dh1) + np.square(gradY[1]/dh2) )
        gradT = np.sqrt( np.square(gradX) + np.square(gradY) )
        
        # Renormalize gradient to yield a signal for "how nonlinear is the transform?"
        (slope_h2x, intercept_h2x), (slope_h1y, intercept_h1y) = self.get_linear_fit()
        slope = np.average((np.abs(slope_h2x), np.abs(slope_h1y))) # nm/mm
        grad_typical = 1/slope # mm/nm
        gradN = gradT/grad_typical
        # gradN = 1 is normal/expected; anything
        # gradN > 1 represents high distortion/compression; h(x) was flat, x(h) is steep
        # gradN < 1 represents high distortion/stretching; h(x) was vertical, x(h) is flat
        gradN = np.where( gradN>=1, gradN, 1/gradN ) - 1
        # Now gradN is [0..+inf]; gradN~1 means gradient was twice as large as expectations
        
        
        # Cyclic delta
        H1c, H2c = self.xy_to_hh(X, Y)
        delta = np.sqrt( np.square(H1-H1c) + np.square(H2-H2c) )
        span = np.average( [ self.h1_range[1]-self.h1_range[0], self.h2_range[1]-self.h2_range[0] ] )
        delta *= 100/span # I.e. we consider a deviation of 1/100 of the range to be 'large'


        # Difference from linear
        XY = np.column_stack( (X.ravel(), Y.ravel()) ) # XY.shape = (N, 2) 
        linear_h1_model, linear_h2_model = self.get_linear_model()
        predict_H1 = linear_h1_model.predict(XY).reshape(H1.shape)
        predict_H2 = linear_h2_model.predict(XY).reshape(H2.shape)
        diff_H1 = (H1-predict_H1)/span
        diff_H2 = (H2-predict_H2)/span

        diff = np.sqrt( np.square(diff_H1) + np.square(diff_H2) )*2
        
        
        # Sum together the different estimators of distortion
        D = np.power(gradN + delta + diff, n_power)
        
        return D




    def plot_distortion(self, name, num=50, plot_buffers=[0.1,0.01,0.11,0.05]):
        import matplotlib as mpl
        mpl.rcParams['mathtext.fontset'] = 'cm'
        import pylab as plt
        
        import time


        self.fig = plt.figure( figsize=(6,6), facecolor='white' )
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )

        self.ax.set_xlabel('$h_{\mathrm{L2}} \, (\mathrm{nm})$', size=14)
        self.ax.set_ylabel('$h_{\mathrm{L1}} \, (\mathrm{nm})$', size=14)
        cmap = mpl.cm.viridis
        
        
        h1_axis = np.linspace(self.h1_range[0], self.h1_range[1], num=num, endpoint=True)
        h2_axis = np.linspace(self.h2_range[0], self.h2_range[1], num=num, endpoint=True)
        H2, H1 = np.meshgrid(h2_axis, h1_axis)
        # Grids are organized as Z[h1_idx, h2_idx]
        
        
        start_time = time.time()
        D = self.distortion_2Darray_map( H1, H2, h1_axis, h2_axis, n_power=1)
        #D = self.distortion_2Darray(H1, H2, h1_axis, h2_axis)
        print('Computing distortion map took {:.1f} s'.format(time.time()-start_time))

        self.im = self.ax.pcolormesh( h2_axis, h1_axis, D, cmap=cmap, shading='flat', vmin=0, vmax=1)
        #, vmin=0, vmax=1)
        self._plot_sample_grid()
        
        #self.ax.axis([xi, xf, yi, yf])
        self.ax.set_aspect('equal')
        self.ax.axis([self.h2_range[0], self.h2_range[1], self.h1_range[0], self.h1_range[1]],)
        
        self.fig.savefig('{}-distortion.png'.format(name), dpi=200)
        plt.close()


    def _plot_sample_grid(self, gnum=4, num=1000, **kwargs):
        conditions = [
            ['k', 2.0, gnum] , # Vertical major lines
            ['k', 0.5, gnum*4] , # Vertical minor lines
            ]
        for c, t, cnum in conditions:
            for x in np.linspace(self.x_range[0], self.x_range[1], num=cnum, endpoint=True):
                y_plot = np.linspace(self.y_range[0], (self.y_range[1]+self.y_range[0])/2, num=num, endpoint=True)
                x_plot = np.ones_like(y_plot)*x
                #self.ax.plot(x_plot, y_plot, '-', color=c, linewidth=t, dashes=[2,1])
                h1, h2 = self.xy_to_hh(x_plot, y_plot)
                self.ax.plot(h2, h1, '-', color=c, linewidth=t, dashes=[2,1])

                y_plot = np.linspace((self.y_range[1]+self.y_range[0])/2, self.y_range[1], num=num, endpoint=True)
                x_plot = np.ones_like(y_plot)*x
                #self.ax.plot(x_plot, y_plot, '-', color=c, linewidth=t)
                h1, h2 = self.xy_to_hh(x_plot, y_plot)
                self.ax.plot(h2, h1, '-', color=c, linewidth=t)


        conditions = [
            ['0.5', 2.0, gnum] , # Horizontal major lines
            ['0.5', 0.5, gnum*4] , # Horizontal minor lines
            ]
        for c, t, cnum in conditions:
            for y in np.linspace(self.y_range[0], self.y_range[1], num=cnum, endpoint=True):
                x_plot = np.linspace(self.x_range[0], (self.x_range[1]+self.x_range[0])/2, num=num, endpoint=True)
                y_plot = np.ones_like(x_plot)*y
                #self.ax.plot(x_plot, y_plot, '-', color=c, linewidth=t)
                h1, h2 = self.xy_to_hh(x_plot, y_plot)
                self.ax.plot(h2, h1, '-', color=c, linewidth=t)

                x_plot = np.linspace((self.x_range[1]+self.x_range[0])/2, self.x_range[1], num=num, endpoint=True)
                y_plot = np.ones_like(x_plot)*y
                #self.ax.plot(x_plot, y_plot, '-', color=c, linewidth=t, dashes=[2,1])
                h1, h2 = self.xy_to_hh(x_plot, y_plot)
                self.ax.plot(h2, h1, '-', color=c, linewidth=t, dashes=[2,1])        


    def plot_trends(self, name='sample', num=100, d_extra=20):
        # Plot the individual conversion equations
        import matplotlib as mpl
        mpl.rcParams['mathtext.fontset'] = 'cm'
        import pylab as plt
        
        self.fig = plt.figure(figsize=(8,5))
        plot_buffers = [0.1,0.01,0.11,0.05]
        gap = 0.05
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 0.5-right_buf-left_buf-gap
        fig_height = 1.0-top_buf-bottom_buf
        self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )
        self.ax2 = self.fig.add_axes( [2*left_buf+fig_width+gap, bottom_buf, fig_width, fig_height] )

        hmin, hmax = None, None
        for y_relative in [0, 0.5, 1.0]:
            # The valid range
            x = np.linspace(self.x_range[0], self.x_range[1], num=num, endpoint=True)
            y = np.ones_like(x)*( self.y_range[0] + (self.y_range[1]-self.y_range[0])*y_relative )
            h1, h2 = self.xy_to_hh(x, y)
            if hmin is None or np.min(h2)<hmin:
                hmin = np.min(h2)
            if hmax is None or np.max(h2)>hmax:
                hmax = np.max(h2)
            
            self.ax.plot(x, h2, '-', color='k', linewidth=2)
            self.ax.axvspan(self.x_range[0], self.x_range[1], color='b', alpha=0.1)
            
            idx = -1
            self.ax.plot([x[idx]], [h2[idx]], 'o', color='k', markersize=5)
            s = 'xy({:.1f}, {:.1f})\nhh[{:.1f}, {:.1f}]'.format(x[idx], y[idx], h1[idx], h2[idx])
            self.ax.text(x[idx], h2[idx], s, size=6, verticalalignment='bottom', horizontalalignment='left')

            idx = 0
            self.ax.plot([x[idx]], [h2[idx]], 'o', color='k', markersize=5)
            s = 'xy({:.1f}, {:.1f})\nhh[{:.1f}, {:.1f}]'.format(x[idx], y[idx], h1[idx], h2[idx])
            self.ax.text(x[idx], h2[idx], s, size=6, verticalalignment='top', horizontalalignment='right')

            # Extended range
            x = np.linspace(self.x_range[0]-d_extra, self.x_range[1]+d_extra, num=num*2, endpoint=True)
            y = np.ones_like(x)*( self.y_range[0] + (self.y_range[1]-self.y_range[0])*y_relative )
            h1, h2 = self.xy_to_hh(x, y)
            self.ax.plot(x, h2, '-', color='k', linewidth=0.5, dashes=[3,3])
        

        self.ax.axhline(hmin, color='b', linewidth=1, dashes=[4,4])
        self.ax.axhline(hmax, color='b', linewidth=1, dashes=[4,4])
        xi, xf, yi, yf = self.ax.axis()
        for h in [hmin, hmax]:
            self.ax.text(xf, h, '{:.1f}'.format(h), size=10, color='b', horizontalalignment='left', verticalalignment='center')

        self.ax.set_xlabel('$x \, (\mathrm{mm}$)', size=14)
        self.ax.set_ylabel('$h_{\mathrm{L2}} \, (\mathrm{nm})$', size=14)
        
        
        hmin, hmax = None, None
        for x_relative in [0, 0.5, 1.0]:
            # The valid range
            y = np.linspace(self.y_range[0], self.y_range[1], num=num, endpoint=True)
            x = np.ones_like(y)*( self.x_range[0] + (self.x_range[1]-self.x_range[0])*x_relative )
            h1, h2 = self.xy_to_hh(x, y)
            if hmin is None or np.min(h1)<hmin:
                hmin = np.min(h1)
            if hmax is None or np.max(h1)>hmax:
                hmax = np.max(h1)
            
            self.ax2.plot(y, h1, '-', color='k', linewidth=2)
            self.ax2.axvspan(self.y_range[0], self.y_range[1], color='b', alpha=0.1)
            
            idx = -1
            self.ax2.plot([y[idx]], [h1[idx]], 'o', color='k', markersize=5)
            s = 'xy({:.1f}, {:.1f})\nhh[{:.1f}, {:.1f}]'.format(x[idx], y[idx], h1[idx], h2[idx])
            self.ax2.text(y[idx], h1[idx], s, size=6, verticalalignment='bottom', horizontalalignment='left')

            idx = 0
            self.ax2.plot([y[idx]], [h1[idx]], 'o', color='k', markersize=5)
            s = 'xy({:.1f}, {:.1f})\nhh[{:.1f}, {:.1f}]'.format(x[idx], y[idx], h1[idx], h2[idx])
            self.ax2.text(y[idx], h1[idx], s, size=6, verticalalignment='top', horizontalalignment='right')

            # Extended range
            y = np.linspace(self.y_range[0]-d_extra, self.y_range[1]+d_extra, num=num*2, endpoint=True)
            x = np.ones_like(y)*( self.x_range[0] + (self.x_range[1]-self.x_range[0])*x_relative )
            h1, h2 = self.xy_to_hh(x, y)
            self.ax2.plot(y, h1, '-', color='k', linewidth=0.5, dashes=[3,3])
        

        self.ax2.axhline(hmin, color='b', linewidth=1, dashes=[4,4])
        self.ax2.axhline(hmax, color='b', linewidth=1, dashes=[4,4])
        xi, xf, yi, yf = self.ax2.axis()
        for h in [hmin, hmax]:
            self.ax2.text(xf, h, '{:.1f}'.format(h), size=10, color='b', horizontalalignment='left', verticalalignment='center')

        self.ax2.set_xlabel('$y \, (\mathrm{mm}$)', size=14)
        self.ax2.set_ylabel('$h_{\mathrm{L1}} \, (\mathrm{nm})$', size=14)        
        
        
        self.fig.savefig('{}-trends.png'.format(name), dpi=200)
        plt.close()
  





# { 'anneal_time': 8, 'pre_anneal_time': 0 }
#para1 = -2.69484324e-03, -7.62126902e-01, -8.25230180e-02,  4.36702251e+01
#para2 = -1.57019486e-02,  2.13823384e-01, -4.42929671e-02,  8.08603651e+01

# { 'anneal_time': 300, 'pre_anneal_time': 15 }
#para1 = -8.04153342e-03, -1.16172684e+00, -9.03572823e-02,  4.71113308e+01
#para2 = -5.69804793e-03,  7.43189749e-01, -3.96324888e-01,  1.10204179e+02

# { 'anneal_time': 30, 'pre_anneal_time': 15 }
#para1 = -9.26422895e-03, -1.33698298e+00, -1.63724597e-01,  5.04172327e+01
#para2 = -8.60228007e-03,  5.53038695e-01, -4.38947684e-01,  1.08959805e+02

# { 'anneal_time': 300, 'pre_anneal_time': 5 }
#para1 = -7.57678073e-03, -1.05516196e+00, -8.66952175e-02,  4.50252003e+01
#para2 = -1.07839722e-02,  2.74579458e-01, -7.46268383e-02,  9.24796652e+01

#sample 15
# { 'anneal_time': 30, 'pre_anneal_time': 5 }
#para1 = -5.11514836e-03, -1.00183291e+00, -1.20151520e-01,  4.89183421e+01
#para2 = -6.66068623e-03,  2.73440303e-01, -1.91366886e-01,  8.22806603e+01


#name = 'sample13'
# { 'anneal_time': 15, 'pre_anneal_time': 5 }
#para1 = -6.83636299e-03, -1.08523640e+00, -2.58277384e-01,  5.94620651e+01
#para2 = -1.86804952e-02, -1.62818801e-01, -2.89007114e-01,  9.28036548e+01

#name = 'sample16'
# { 'anneal_time': 1200, 'pre_anneal_time': 5 }
#para1 = 3.45665596e-03, -9.48340775e-01, -2.15781785e-01,  5.16824292e+01
#para2 = 1.34131761e-03,  8.57767669e-01, -1.31882004e-01,  9.43731748e+01

#name = 'sample7'
# { 'anneal_time': 1200, 'pre_anneal_time': 15 }
#para1 = -1.09572715e-02, -1.33776535e+00, -1.55529610e-01,  4.88982389e+01
#para2 = -5.04082592e-03,  5.38322523e-01, -4.31686426e-01,  1.01201989e+02
#x bounds: -2.0 to 52.0
#y bounds: -2.0 to 52.0
#h1 bounds: 38.1 to 89.1
#h2 bounds: 37.1 to 103.1

#name = 'sample6'
# { 'anneal_time': 15, 'pre_anneal_time': 0 }
#para1 = -5.47709055e-03, -1.11144308e+00, -9.34935198e-02,  4.44242154e+01
#para2 = -1.94166224e-02, -7.92471986e-02, -1.80287725e-01,  8.95367666e+01
#x bounds: -2.0 to 52.0
#y bounds: -2.0 to 52.0
#h1 bounds: 37.3 to 87.6
#h2 bounds: 31.8 to 90.0

name = 'sample11'
# { 'anneal_time': 0, 'pre_anneal_time': 15 }
para1 = -6.66921517e-03, -1.10472825e+00, -2.75380237e-01,  5.74740059e+01
para2 = -2.07761774e-02, -2.35510686e-01, -1.91280054e-01,  8.38146998e+01
#x bounds: -2.0 to 52.0
#y bounds: -2.0 to 52.0
#h1 bounds: 40.9 to 97.4
#h2 bounds: 29.9 to 84.8

#name = 'sample12'
# { 'anneal_time': 0, 'pre_anneal_time': 5 }
#para1 = -4.92994429e-03, -9.66557810e-01, -2.35087869e-01,  5.07872174e+01
#para2 = -7.93071791e-03,  5.16553479e-01, -1.17074592e-01,  9.25285055e+01
#x bounds: -2.0 to 52.0
#y bounds: -2.0 to 52.0
#h1 bounds: 36.6 to 88.2
#h2 bounds: 38.1 to 93.8





#name = 'LonC_sample2' 
# { 'anneal_time': 0, 'pre_anneal_time': 0 }
#para1 = 1.00790506e-03, -5.53822740e-01, -6.28534755e-02,  4.11339095e+01
#para2 = -1.47514140e-02,  5.28820379e-01,  1.57572692e-01,  7.91861305e+01
#x bounds: -2.0 to 52.0
#y bounds: -2.0 to 52.0
#h1 bounds: 36.8 to 72.8
#h2 bounds: 11.5 to 88.4


#name = 'LonC_sample3'
# { 'anneal_time': 15, 'pre_anneal_time': 0 }
#para1 = -3.32167067e-03, -9.26626393e-01, -2.16140860e-01,  5.47694514e+01
#para2 = -2.10060483e-02,  4.49004236e-01,  1.78344541e-01,  9.32566125e+01
#x bounds: -2.0 to 52.0
#y bounds: -2.0 to 52.0
#h1 bounds: 41.7 to 94.4
#h2 bounds: 12.8 to 103.3

#name = 'LonC_sample1'
# { 'anneal_time': 8, 'pre_anneal_time': 0 }
#para1 = -6.27355605e-03, -1.03319924e+00, -1.38196229e-01,  5.10194140e+01
#para2 = -2.73036983e-02,  2.90190044e-01,  2.46008824e-01,  8.94386732e+01
#x bounds: -2.0 to 52.0
#y bounds: -2.0 to 52.0
#h1 bounds: 41.7 to 88.1
#h2 bounds: 0.0 to 102.7


#name = 'sample5'
# { 'anneal_time': 0, 'pre_anneal_time': 0 }
#para1 = -1.11096877e-03, -8.02799083e-01, -7.98218309e-02,  4.34740265e+01
#para2 = -1.98990675e-02, -2.40318834e-01, -5.27970432e-02,  7.11330200e+01
#x bounds: -2.0 to 52.0
#y bounds: -2.0 to 52.0
#h1 bounds: 37.7 to 82.4
#h2 bounds: 27.1 to 72.0



coord = CoordinateTransformation(para1, para2, x_range=[-2, 52], y_range=[-2, 52], method='fsolve')
#coord.hh_to_xy_linear_model(70, 80)

#print(coord.hh_to_xy(57.79,93.0))
#coord.confirm_working()
#coord.determine_bounds()
#coord.plot_trends(name)
#coord.plot_maps(name)

#coord.plot_distortion(name)



