#!/usr/bin/python

import os
import time, datetime
from databroker import list_configs
#print( list_configs() )
from databroker import Broker





# Select the experiment to consider
beamline = 'cms'
start_time, stop_time = '2019-01-01', '2019-01-30'
#alias_dir = '/GPFS/xf11bm/data/2019_1/UShell'
alias_dir = os.path.abspath('./')



# Select what kinds of scans to consider
measure_type = 'measure'

# Select the meta-data (md) to output
md_required = ['scan_id']
#md_required = ['scan_id', 'sample_x', 'sample_y', 'sample_th']
#md_required = ['scan_id', 'motor_SAXSx', 'motor_SAXSy', 'motor_DETx', 'motor_WAXSx']

md_optional = []
#md_optional = ['sample_clock', 'sample_temperature']
#md_optional = ['mfc', 'film_thickness']
#md_optional = ['sample_x', 'sample_y', 'sample_th', 'sample_motor_x', 'sample_motor_y', 'sample_motor_th', 'sample_clock', 'T_actual']







# Additional settings
verbosity = 3
outfile = 'manifest.txt'
delimeter = '\t'
missing_md = '-'




if verbosity>=1:
    print('\n========================================')
    print('Generating manifest...')
    print('========================================')

execution_start = time.time()
db = Broker.named(beamline)

if verbosity>=1:
    print("Scans (of type '{}') in directory:".format(measure_type))
    print("    {}".format(alias_dir))
#headers = db(start_time=start_time, stop_time=stop_time, measure_type=measure_type, experiment_alias_directory=alias_dir)
headers = db(since=start_time, until=stop_time, measure_type=measure_type, experiment_alias_directory=alias_dir)


if verbosity>=3:
    header_list = list(headers)
    num_scans = len(header_list)
    print('{} measurements'.format(num_scans))
else:
    num_scans = None


with open(outfile, 'w') as fout:
    
    fields = ['datetime', 'timestamp', 'name']
    for mdc in md_required:
        fields.append(mdc)
    for mdc in md_optional:
        fields.append(mdc)
    headline = delimeter.join(fields)
    fout.write('#{}\n'.format(headline))
    
    for ih, header in enumerate(headers):
        
        if verbosity>=2:
            if ih%100==0:
                if num_scans is None:
                    print(' Record {}'.format(ih))
                else:
                    print(' Record {} of {} ({:.1f}% done after {:.2f}s)'.format(ih, num_scans, 100.0*ih/num_scans, time.time()-execution_start ))
        
        
        
        if verbosity>=6:
            # Print out diagnostic info
            #print(header)
            #print(header.table()['det'].mean())
            for k, v in header.items():
                print(k, v)

        section = header['start']
        
        if verbosity>=5:
            for k, v in section.items():
                print(k, v)

        unixtime = float( section['time'] )
        time_str = datetime.datetime.fromtimestamp(unixtime).strftime('%Y-%m-%d %H:%M:%S')

        #name = section['sample_savename']
        name = section['filename']
        #seqID = section['detector_sequence_ID']
        #sID = section['sample_measurement_ID']

        fields = [time_str, unixtime, name]

        for mdc in md_required:
            fields.append(section[mdc])

        for mdc in md_optional:
            if mdc in section:
                fields.append(section[mdc])
            else:
                fields.append(missing_md)
                


        fields = [str(f) for f in fields]
        line = delimeter.join(fields)

        if verbosity>=4:
            print(line)
        fout.write('{}\n'.format(line))


if verbosity>=1:
    print('========================================')
    print('Done. Execution took {:.2f} s'.format(time.time()-execution_start))
    print('========================================\n')
