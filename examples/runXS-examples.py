#test

class metadata_extract(tools.Protocol):
    
    def __init__(self, name='metadata_extract', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.dat'
        self.run_args = {}
        self.run_args.update(kwargs)
    
        
    @Protocols.run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        infile = data.infile
        f = tools.Filename(infile)
        filepath, filename, filebase, ext = f.split()
        
        results['infile'] = infile
        results['filepath'] = filepath
        results['filename'] = filename
        results['filebase'] = filebase
        results['fileext'] = ext
        
        results['sample_name'] = data.name
        results['file_ctime'] = os.path.getctime(infile)
        results['file_modification_time'] = os.path.getmtime(infile)
        results['file_access_time'] = os.path.getatime(infile)
        results['file_size'] = os.path.getsize(infile)
        
        
        patterns = [
                    #['theta', '.+_th(\d+\.\d+)_.+'] ,
                    #['annealing_temperature', '.+_T(\d+\.\d\d\d)C_.+'] ,
                    #['annealing_time', '.+_(\d+\.\d)s_T.+'] ,
                    ['x_position', '.+_x(-?\d+\.\d+)_.+'] ,
                    ['y_position', '.+_y(-?\d+\.\d+)_.+'] ,
                    ['exposure_time', '.+_(\d+\.\d+)c_\d+_saxs.+'] ,
                    ['sequence_ID', '.+_(\d+)_saxs.+'] ,
                    ]
                    
        
        for pattern_name, pattern_string in patterns:
            pattern = re.compile(pattern_string)
            m = pattern.match(filename)
            if m:
                if run_args['verbosity']>=5:
                    print('  matched: {} = {}'.format(pattern_name, m.groups()[0]))
                results[pattern_name] = float(m.groups()[0])
        
        
        #outfile = self.get_outfile(data.name, output_dir)
        
        
        return results
    


# Results
################################################################################
class Results(object):
    '''Simple object to help extract result values from a bunch of xml files.'''
    
    def __init__(self):
        
        import xml.etree.ElementTree as etree
        #from lxml import etree
        #import xml.dom.minidom as minidom
        
        self.etree = etree
        
        
    def extract_save_txt(self, outfile, infiles, protocol, result_names):
        
        results = self.extract(infiles, protocol, result_names)
        
        with open(outfile, 'w') as fout:
            
            #fout.write('
            
            header = '#filename\t{}\n'.format('\t'.join(result_names))
            fout.write(header)
            
            for result in results:
                line = '{}\n'.format('\t'.join(str(r) for r in result))
                fout.write(line)
                
    
    def extract(self, infiles, protocol, result_names):
        '''Extract the specified results-values (for the given protocol), from
        the specified files. The most recent run of the protocol is used.'''
        
        results = []
        
        for i, infile in enumerate(infiles):
            
            if i%100==0:
                print( 'Extracting file {} ({:.1f}% done)'.format(i+1, 100.*i/len(infiles)) )
            line = [infile]
            
            result_names_e, results_e = self.extract_results_from_xml(infile, protocol)
            
            if len(result_names_e)>0:
            
                for result_name in result_names:
                    idx = result_names_e.index(result_name)
                    line.append(results_e[idx])
                    
                results.append(line)
                
                
        return results
            
            
    def extract_multi(self, infiles, extractions):
        
        results = []
        for i, infile in enumerate(infiles):
            results.append( [infile] )
        
        for i, infile in enumerate(infiles):
            
            if i%100==0:
                print( 'Extracting file {} ({:.1f}% done)'.format(i+1, 100.*i/len(infiles)) )
            
            for protocol, result_names in extractions:
            
                result_names_e, results_e = self.extract_results_from_xml(infile, protocol)
                
                for result_name in result_names:
                    if result_name in result_names_e:
                        idx = result_names_e.index(result_name)
                        results[i].append(results_e[idx])
                    else:
                        results[i].append('-')
                    
                #results.append(line)
                
        return results            
        
    def extract_results_from_xml(self, infile, protocol):
        
        #parser = self.etree.XMLParser(remove_blank_text=True)
        parser = self.etree.XMLParser()
        root = self.etree.parse(infile, parser).getroot()

        # Get the latest protocol
        element = root
        children = [child for child in element if child.tag=='protocol' and child.get('name')==protocol]
        children_v = [float(child.get('end_timestamp')) for child in element if child.tag=='protocol' and child.get('name')==protocol]

        result_names = []
        results = []
        
        if len(children_v)>0:
            idx = np.argmax(children_v)
            protocol = children[idx]
            
            # In this protocol, get all the results (in order)
            element = protocol
            children = [child for child in element if child.tag=='result']
            children_v = [child.get('name') for child in element if child.tag=='result']
            
            idx = np.argsort(children_v)
            #result_elements = np.asarray(children)[idx]
            result_elements = [children[i] for i in idx]
            
            for element in result_elements:
                
                #print( element.get('name') )
                
                if element.get('value') is not None:
                    result_names.append(element.get('name'))
                    try:
                        results.append(float(element.get('value')))
                    except ValueError:
                        results.append(element.get('value'))
                    
                    if element.get('error') is not None:
                        result_names.append(element.get('name')+'_error')
                        results.append(float(element.get('error')))
                    
                elif element.get('type') is not None and element.get('type')=='list':
                    
                    # Elements of the listextract_save_txt
                    children = [child for child in element if child.tag=='element']
                    children_v = [int(child.get('index')) for child in element if child.tag=='element']
                    #print(children_v)
                    
                    # Sorted
                    idx = np.argsort(children_v)
                    children = [children[i] for i in idx]
                    
                    # Append values
                    for child in children:
                        #print( child.get('index') )
                        result_names.append('{}_{}'.format(element.get('name'), child.get('index')))
                        results.append(float(child.get('value')))
                        
                
                else:
                    print('    Errror: result has no usable data ({})'.format(element))
            
        
        return result_names, results
                
    
    # End class Results(object)
    ########################################
    



class update_autonomous_data(Protocols.Protocol):
    
    def __init__(self, name='autonomous', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.npy'
        self.run_args = {
                        'crop' : None,
                        'shift_crop_up' : 0.0,
                        'blur' : 2.0,
                        'resize' : 0.2,
                        'ztrim' : [0.05, 0.005]
                        }
        self.run_args.update(kwargs)
        

    @tools.run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        import time
        time.sleep(1) # Kludge to avoid corrupting XML file?
        
        
        # Compile results
        results_xml_file = self.get_outfile(data.name, './results/', ext='.xml')
        
        if run_args['verbosity']>=5:
            print('    Extracting results from: {}'.format(results_xml_file))
        
        
        
        infiles = [results_xml_file]
        #extractions = [ [ 'metadata_extract', ['x_position','y_position'] ] ,
                    #['circular_average_q2I_fit', ['fit_peaks_d0', 'fit_peaks_grain_size', 'fit_peaks_prefactor1'] ],
                    #]
        extractions = [ [ 'metadata_extract', ['x_position','y_position'] ] ,
                    ['circular_average_q2I_fit', ['fit_peaks_prefactor1'] ],
                    ]

        
        extracted_results = Results().extract_multi(infiles, extractions)
        data_vector = np.asarray(extracted_results)[0,1:]
        data_vector = data_vector.astype(np.float)


        
        if 'file_extension' in run_args and run_args['file_extension'] is not None:
            outfile = self.get_outfile(data.name, output_dir, ext=run_args['file_extension'])
        else:
            outfile = self.get_outfile(data.name, output_dir)
        
        #print(data.stats())
        
        np.save(outfile, data_vector)
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'data vector (npy)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
        
        
        
        # Update the file for SMART to analyze
        outfile = '../autonomous/data_new/data.npy'
        np.save(outfile, data_vector)
        
        
        return results    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


class circular_average_subtract(Protocols.circular_average):

    def __init__(self, name='circular_average_subtract', **kwargs):

        self.name = self.__class__.__name__ if name is None else name

        self.default_ext = '.png'
        self.run_args = {}
        self.run_args.update(kwargs)

    @tools.run_default
    def run(self, data, output_dir, **run_args):

        results = {}

        if 'dezing' in run_args and run_args['dezing']:
            data.dezinger(sigma=3, tol=100, mode='median', mask=True, fill=False)


        line = data.circular_average_q_bin(error=True)
        #line.smooth(2.0, bins=10)


        background_filename = './circular_average/blank_300.00s_314934_saxs.dat'
        background_data = np.loadtxt(background_filename)
        background = background_data[:,2]

        line.y -= background


        outfile = self.get_outfile(data.name, output_dir)

        try:
            line.plot(save=outfile, **run_args)
        except ValueError:
            pass

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)


        return results




    
class circular_average_q2I_fit(Protocols.circular_average_q2I):

    def __init__(self, name=None, **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {'show_region' : False ,
                         'plot_range' : [None, None, 0, None] ,
                         'auto_plot_range_fit' : True ,
                         }
        self.run_args.update(kwargs)
    
        
    @tools.run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        line = data.circular_average_q_bin(error=True)
        
        #line.y *= np.square(line.x)
        #line.y *= np.power(line.x, 0.75)
        line.y *= np.abs(line.x)
        
        
        line.y_label = 'q^n *I(q)'
        line.y_rlabel = r'$q^n I(q) \, (\mathrm{ \AA^{-n} \times counts/pixel})$'
        
        
        #outfile = self.get_outfile(data.name, output_dir, ext='_q2I{}'.format(self.default_ext))
        #line.plot(save=outfile, show=False, **run_args)
        outfile = self.get_outfile(data.name, output_dir, ext='_q2I.dat')
        line.save_data(outfile)        

        
        # Fit data
        #if 'fit_range' in run_args:
            #line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
            #line.trim(run_args['fit_range'][0], run_args['fit_range'][1])
        if 'trim_range' in run_args:
            line.trim(run_args['trim_range'][0], run_args['trim_range'][1])
        
        lm_result, fit_line, fit_line_extended = self._fit_peaks(line, **run_args)
        
        
        # Save fit results
        fit_name = 'fit_peaks'
        prefactor_total = 0
        for param_name, param in lm_result.params.items():
            results['{}_{}'.format(fit_name, param_name)] = { 'value': param.value, 'error': param.stderr, }
            if 'prefactor' in param_name:
                prefactor_total += np.abs(param.value)
            
        results['{}_prefactor_total'.format(fit_name)] = prefactor_total
        results['{}_chi_squared'.format(fit_name)] = lm_result.chisqr/lm_result.nfree
        
        # Calculate some additional things
        d = 0.1*2.*np.pi/results['{}_x_center1'.format(fit_name)]['value']
        results['{}_d0'.format(fit_name)] = d
        xi = 0.1*(2.*np.pi/np.sqrt(2.*np.pi))/results['{}_sigma1'.format(fit_name)]['value']
        results['{}_grain_size'.format(fit_name)] = xi       

        if False:
            d = 0.1*2.*np.pi/results['{}_x_center2'.format(fit_name)]['value']
            results['{}_d02'.format(fit_name)] = d
            xi = 0.1*(2.*np.pi/np.sqrt(2.*np.pi))/results['{}_sigma2'.format(fit_name)]['value']
            results['{}_grain_size2'.format(fit_name)] = xi              
        
        # Plot and save data
        class DataLines_current(DataLines):
            
            def _plot_extra(self, **plot_args):
                
                xi, xf, yi, yf = self.ax.axis()
                v_spacing = (yf-yi)*0.08
                
                yp = yf
                s = '$a = \, {:.3f}$'.format(self.results['fit_peaks_prefactor1']['value'])
                self.ax.text(xi, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='left')
                
                
                yp -= v_spacing
                s = '$q_0 = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(self.results['fit_peaks_x_center1']['value'])
                self.ax.text(xi, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='left')

                yp -= v_spacing
                s = r'$d_0 \approx \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_d0'])
                self.ax.text(xi, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='left')

                yp -= v_spacing
                s = '$\sigma = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(self.results['fit_peaks_sigma1']['value'])
                self.ax.text(xi, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='left')
                
                yp -= v_spacing
                s = r'$\xi \approx \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_grain_size'])
                self.ax.text(xi, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='left')


                if False:
                    yp = yf
                    s = '$a = \, {:.3f}$'.format(self.results['fit_peaks_prefactor2']['value'])
                    self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')
                    
                    yp -= v_spacing
                    s = '$q_0 = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(self.results['fit_peaks_x_center2']['value'])
                    self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')

                    yp -= v_spacing
                    s = r'$d_0 \approx \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_d02'])
                    self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')

                    yp -= v_spacing
                    s = '$\sigma = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(self.results['fit_peaks_sigma2']['value'])
                    self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')
                    
                    yp -= v_spacing
                    s = r'$\xi \approx \, {:.1f} \, \mathrm{{nm}}$'.format(self.results['fit_peaks_grain_size2'])
                    self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')


        
        lines = DataLines_current([line, fit_line, fit_line_extended])
        lines.copy_labels(line)
        lines.results = results

        outfile = self.get_outfile(data.name+'-fit', output_dir, ext='.png')
        
        # Tweak the plotting range for the fit-plot
        run_args_cur = run_args.copy()
        if run_args['auto_plot_range_fit']:
            run_args_cur['plot_range'] = [ run_args['plot_range'][0] , run_args['plot_range'][1] , run_args['plot_range'][2] , run_args['plot_range'][3] ]
            if 'fit_range' in run_args_cur:
                span = abs(run_args['fit_range'][1]-run_args_cur['fit_range'][0])
                run_args_cur['plot_range'][0] = run_args['fit_range'][0]-span*0.25
                run_args_cur['plot_range'][1] = run_args_cur['fit_range'][1]+span*0.25
            
            run_args_cur['plot_range'][2] = 0
            run_args_cur['plot_range'][3] = max(fit_line.y)*1.3
        
        
        try:
            #lines.plot(save=outfile, error_band=False, ecolor='0.75', capsize=2, elinewidth=1, **run_args)
            lines.plot(save=outfile, **run_args_cur)
        except ValueError:
            pass


        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.x_err = None
        line.y_err = None
        line.save_data(outfile)
        
        return results
                       

    def _fit_peaks(self, line, num_curves=1, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self.fit_peaks(line, **run_args)

        line_full = line
        if 'fit_range' in run_args:
            line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
        import lmfit
        
        def model(v, x):
            
            # Linear background
            m = v['m']*x + v['b']
            # Power-law background
            m += v['qp']*np.power( np.abs(x), v['qalpha'] )
            
            # Gaussian peaks
            for i in range(num_curves):
                m += v['prefactor{:d}'.format(i+1)]*np.exp( -np.square(x-v['x_center{:d}'.format(i+1)])/(2*(v['sigma{:d}'.format(i+1)]**2)) )
            
            return m
        
        def func2minimize(params, x, data):
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()

        m = (line.y[-1]-line.y[0])/(line.x[-1]-line.x[0])
        b = line.y[0] - m*line.x[0]

        xs = np.abs(line.x)
        ys = line.y
        qalpha = (np.log(ys[0])-np.log(ys[-1]))/(np.log(xs[0])-np.log(xs[-1]))
        qp = np.exp( np.log(ys[0]) - qalpha*np.log(xs[0]) )

        if True:
            # Linear background
            params.add('m', value=m, min=0, max=abs(m)*+4, vary=False)
            params.add('b', value=b, min=np.max(line.y)*-100, max=np.max(line.y)*100, vary=False)
            
            params.add('qp', value=0, vary=False)
            params.add('qalpha', value=1.0, vary=False)
            
        else:
            # Power-law background
            params.add('m', value=0, vary=False)
            params.add('b', value=0, vary=False)
            
            params.add('qp', value=qp, vary=False)
            params.add('qalpha', value=qalpha, vary=False)
            
        
        xspan = np.max(line.x) - np.min(line.x)
        xpeak, ypeak = line.target_y(np.max(line.y))
        
        # Best guess for peak position
        if True:
            # Account for power-law scaling (Kratky-like)
            xs = np.asarray(line.x)
            ys = np.asarray(line.y)
            
            ys = ys*np.power( np.abs(xs), np.abs(qalpha) ) # Kratky-like
            
            # Sort
            indices = np.argsort(ys)
            x_sorted = xs[indices]
            y_sorted = ys[indices]
            
            target = np.max(ys)

            # Search through y for the target
            idx = np.where( y_sorted>=target )[0][0]
            xpeak = x_sorted[idx]
            ypeak = y_sorted[idx]
            
            xpeak, ypeak = line.target_x(xpeak)
                                 

        prefactor = ypeak - ( m*xpeak + b )
        sigma = 0.05*xspan
        
        for i in range(num_curves):
            params.add('prefactor{:d}'.format(i+1), value=prefactor, min=0, max=np.max(line.y)*1.5, vary=False)
            params.add('x_center{:d}'.format(i+1), value=xpeak, min=np.min(line.x), max=np.max(line.x), vary=False)
            params.add('sigma{:d}'.format(i+1), value=sigma, min=0, max=xspan*0.75, vary=False)


        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if False:
            # Fit only the peak positions
            params['x_center1'].vary = True
            lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
            params['x_center2'].vary = True
            lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        
        if False:
            # Tweak peak position
            lm_result.params['sigma1'].vary = False
            lm_result.params['x_center1'].vary = True
            lm_result = lmfit.minimize(func2minimize, lm_result.params, args=(line.x, line.y))

        if False:
            # Background
            lm_result.params['m'].vary = True
            lm_result.params['b'].vary = True
            lm_result = lmfit.minimize(func2minimize, lm_result.params, args=(line.x, line.y))
        
        #lmfit.report_fit(lm_result.params)
        #print(lm_result.chisqr)
        
        if True:
            # Relax entire fit
            lm_result.params['m'].vary = True
            lm_result.params['b'].vary = True
            #lm_result.params['qp'].vary = True
            #lm_result.params['qalpha'].vary = True

            for i in range(num_curves):
                #lm_result.params['prefactor{:d}'.format(i+1)].value = lm_result.params['prefactor{:d}'.format(i+1)].value*1.0001 + 1e-14
                
                lm_result.params['prefactor{:d}'.format(i+1)].vary = True
                lm_result.params['sigma{:d}'.format(i+1)].vary = True
                lm_result.params['x_center{:d}'.format(i+1)].vary = True
            
            lm_result = lmfit.minimize(func2minimize, lm_result.params, args=(line.x, line.y))
            lm_result = lmfit.minimize(func2minimize, lm_result.params, args=(line.x, line.y), method='nelder')
            
            
            
        #lmfit.report_fit(lm_result.params)
        #print(lm_result.chisqr)
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'marker':None, 'linewidth':4.0})
        
        #fit_x = np.linspace(np.min(line_full.x), np.max(line_full.x), num=200)
        #fit_x = np.linspace(np.average( [np.min(line_full.x), np.min(line.x)] ), np.average( [0, np.max(line.x)] ), num=200)
        fit_x = np.linspace(np.max(line.x)*0.1, np.max(line.x)*1.2, num=200)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'b', 'alpha':0.5, 'marker':None, 'linewidth':2.0})        

        return lm_result, fit_line, fit_line_extended
    

    
    
    

    
class update_autonomous_data(Protocols.Protocol):
    
    def __init__(self, name='autonomous', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.npy'
        self.run_args = {
                        'crop' : None,
                        'shift_crop_up' : 0.0,
                        'blur' : 2.0,
                        'resize' : 0.2,
                        'ztrim' : [0.05, 0.005]
                        }
        self.run_args.update(kwargs)
        

    @tools.run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        import time
        time.sleep(1) # Kludge to avoid corrupting XML file?
        
        
        # Compile results
        results_xml_file = self.get_outfile(data.name, './results/', ext='.xml')
        
        if run_args['verbosity']>=5:
            print('    Extracting results from: {}'.format(results_xml_file))
        
        
        
        infiles = [results_xml_file]
        #extractions = [ [ 'metadata_extract', ['x_position','y_position'] ] ,
                    #['circular_average_q2I_fit', ['fit_peaks_d0', 'fit_peaks_grain_size', 'fit_peaks_prefactor1'] ],
                    #]
        extractions = [ [ 'metadata_extract', ['x_position','y_position'] ] ,
                    ['circular_average_q2I_fit', ['population_0', 'population_1', 'population_2'] ],
                    ]

        
        extracted_results = Results().extract_multi(infiles, extractions)
        data_vector = np.asarray(extracted_results)[0,1:]
        data_vector = data_vector.astype(np.float)


        
        if 'file_extension' in run_args and run_args['file_extension'] is not None:
            outfile = self.get_outfile(data.name, output_dir, ext=run_args['file_extension'])
        else:
            outfile = self.get_outfile(data.name, output_dir)
        
        #print(data.stats())
        
        np.save(outfile, data_vector)
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'data vector (npy)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
        
        
        
        # Update the file for SMART to analyze
        #outfile = '../autonomous/data_new/data.npy'
        outfile = '/GPFS/xf11bm/data/2018_2/MNoack_codes/June2018/smart1/smart/G-SMART/data/new_experiment_result/experiment_result.npy'
        
        
        np.save(outfile, data_vector)
        
        #if os.path.isfile(outfile):
            ## File exists; append data instead of over-writing
            #old_data_vector = np.load(outfile)
            #data_vector = np.append( old_data_vector, [data_vector], axis=0 )
            #np.save(outfile, data_vector)
            
        #else:
            #np.save(outfile, [data_vector])
        
        
        return results
