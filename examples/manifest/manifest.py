#!/usr/bin/python

import time, datetime
from databroker import list_configs
#print( list_configs() )
from databroker import Broker





# Select the experiment to consider
beamline = 'cms'
start_time, stop_time = '2017-11-27', '2017-12-01'
user = 'RVerduzco'
alias_dir = '/GPFS/xf11bm/data/2017_3/{}'.format(user)

# Select what kinds of scans to consider
measure_type = 'measure'

# Select the meta-data (md) to output
md_required = []
md_optional = ['mfc', 'film_thickness']




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

print("Scans (of type '{}') in directory:".format(measure_type))
print("    {}".format(alias_dir))
headers = db(start_time=start_time, stop_time=stop_time, measure_type=measure_type, experiment_alias_directory=alias_dir)


if verbosity>=3:
    header_list = list(headers)
    num_scans = len(header_list)
    print('{} measurements'.format(num_scans))
else:
    num_scans = None


with open(outfile, 'w') as fout:
    
    fields = ['datetime', 'name']
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
        
        #print(header.table()['det'].mean())

        section = header['start']

        unixtime = float( section['time'] )
        time_str = datetime.datetime.fromtimestamp(unixtime).strftime('%Y-%m-%d %H:%M:%S')

        #name = section['sample_savename']
        name = section['filename']
        #seqID = section['detector_sequence_ID']
        #sID = section['sample_measurement_ID']

        fields = [time_str, name]

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
