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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
class circular_average_sum(Protocols.circular_average):

    def __init__(self, name='circular_average_sum', **kwargs):

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

        line_sub = line.sub_range(run_args['sum_range'][0], run_args['sum_range'][1])
        results['values_sum'] = np.sum(line_sub.y)

        outfile = self.get_outfile(data.name, output_dir)


        # Plot and save data
        class DataLines_current(DataLines):

            def _plot_extra(self, **plot_args):

                xi, xf, yi, yf = self.ax.axis()
                v_spacing = (yf-yi)*0.10

                yp = yf
                s = '$S = \, {:.1f} $'.format(self.results['values_sum'])
                self.ax.text(xf, yp, s, size=20, color='b', verticalalignment='top', horizontalalignment='right')




        lines = DataLines_current([line])
        lines.copy_labels(line)
        lines.results = results

        try:
            lines.plot(save=outfile, **run_args)
        except ValueError:
            pass

        outfile = self.get_outfile(data.name, output_dir, ext='.dat')
        line.save_data(outfile)

        # TODO: Fit 1D data

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
        
        line.y *= np.square(line.x)
        line.y_label = 'q^2*I(q)'
        line.y_rlabel = '$q^2 I(q) \, (\AA^{-2} \mathrm{counts/pixel})$'
        
        
        outfile = self.get_outfile(data.name, output_dir, ext='_q2I{}'.format(self.default_ext))
        line.plot(save=outfile, show=False, **run_args)
        
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
        
        if self._peak_snr>0.05:
            # Calculate some additional things
            d = 0.1*2.*np.pi/results['{}_x_center1'.format(fit_name)]['value']
            results['{}_d0'.format(fit_name)] = d
            xi = 0.1*(2.*np.pi/np.sqrt(2.*np.pi))/results['{}_sigma1'.format(fit_name)]['value']
            results['{}_grain_size'.format(fit_name)] = xi       
            
        else:
            results['{}_d0'.format(fit_name)] = 0
            results['{}_grain_size'.format(fit_name)] = 0
        
        
        # Plot and save data
        class DataLines_current(DataLines):
            
            def _plot_extra(self, **plot_args):
                
                xi, xf, yi, yf = self.ax.axis()
                v_spacing = (yf-yi)*0.10
                
                yp = yf
                s = '$q_0 = \, {:.4f} \, \mathrm{{\AA}}^{{-1}}$'.format(self.results['fit_peaks_x_center1']['value'])
                self.ax.text(xf, yp, s, size=20, color='b', verticlass update_autonomous_data(Protocols.Protocol):
    
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
    