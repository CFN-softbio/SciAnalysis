# The measure_queue should have been loaded from CustomQueue.py
# This should be done prior to loading user.py so that you do
# not keep re-creating the queue (which will cause problems
# since the old queue is already connected to the port). Thus
# in some earlier "%run -i" you should have run code like:
'''
# Connect to ZeroMQ
queue_PATH='../'
#queue_PATH='/nsls2/xf11bm/data/2020_3/MNoack/'
queue_PATH in sys.path or sys.path.append(queue_PATH)
from CustomQueue import *
measure_queue = Queue_measure()
'''

# Note that the use of a relative path, and the fact that
# CustomQueue saves some local fiels, means that you need
# to change to the desired user folder *before* running
# your BlueSky session.



# Example class using autonomous method
# Search for "# TOCHANGE" below for places where beamline-specific
# and experiment-specific assumptions are being made.

class SampleTSAXS(SampleTSAXS_Generic):
    
    def measureAutonomous(self,  exposure_time=0.5, extra=None, max_loop=None, max_measurements=None, prefix='measureAutonomous > ', verbosity=3, **md):
        '''Measure points in a loop, relying on an external queue to specify what
        position to actually measure. If the 'position' is not (x,y) sample coordinates,
        then you will have to add code to do the appropriate coordinate conversion,
        or trigger the right beamline motors/components.'''
        
        iloop = 0
        while True: # Loop forever
            iloop +=1
            if (max_loop is not None) and iloop>max_loop:
                if verbosity>=2:
                    print('{}Reached target maximum iterations ({}). Done.'.format(prefix, max_loop))
                return        
        
            if verbosity>=3:
                print('{}Waiting for AE command on queue...'.format(prefix))
                
            commands = measure_queue.get() # Get measurement command from queue
            num_to_measure = int( sum( [1.0 for command in commands if command['measured'] is False] ) )

            if (max_measurements is not None) and len(commands)>max_measurements:
                if verbosity>=2:
                    print('{}Reached target maximum measurements ({}). Done.'.format(prefix, max_measurements))
                return

            if verbosity>=3:
                print('{}Received command to measure {} points'.format(prefix, num_to_measure))

            
            imeasure = 0
            for icommand, command in enumerate(commands):
                if verbosity>=5:
                    print('{}Considering point {}/{}'.format(prefix, icommand, len(commands)))
                    
                if not command['measured']:
                    imeasure += 1
                    if verbosity>=3:
                        print('{}Measuring point {}/{:d}'.format(prefix, imeasure, num_to_measure))

                    start_time = time.time()
                    
                    
                    
                    ########################################
                    # Move to point
                    ########################################
                    # Here you should define the beamline changes needed to go 
                    # to the desired position. (You shouldn't in general need 
                    # to change code outside of this block.)
                    ########################################
                    
                    # TOCHANGE
                    #a = 25.0 # mm/s^2 (accelaration during the printing)
                    #print_speed = command['position']['print_speed']
                    #x_pos = np.square( np.abs(print_speed) )/(2*a)
                    x_pos = command['position']['x_position']
                    y_pos = command['position']['y_position']

                    if verbosity>=3:
                        print('{}Driving to point {}/{}; (x,y) = ({:.3f}, {:.3f})'.format(prefix, imeasure, num_to_measure, x_pos, y_pos))

                    self.xabs(x_pos)
                    self.yabs(y_pos)

                    # Wait for motors to finish moving
                    #while piezo.x.moving==True:
                        #time.sleep(0.2)
                    #while piezo.y.moving==True:
                        #time.sleep(0.2)
                        
                    
                    command['position']['x_position'] = self.xpos(verbosity=0)
                    command['position']['y_position'] = self.ypos(verbosity=0)
                    #command['position']['print_speed'] = np.sqrt( max(0.0, 2*a*self.xpos(verbosity=0)) )

                    
                    ########################################
                    
                    
                    self.measure(exposure_time=exposure_time, extra=extra, **md)
                    cost_time = time.time() - start_time 
                    
                    header = db[-1] # The most recent measurement
                    command['filename'] = '{}'.format(header.start['filename'][:-1]) # TOCHANGE
                    #command['filename'] = '{}'.format(header.start['filename']) # TOCHANGE
                    command['uid'] = '{}'.format(header.start['uid'])
                    command['cost'] = cost_time
                    command['x_position'] = self.xpos(verbosity=0)
                    command['y_position'] = self.ypos(verbosity=0)
                    command['measured'] = True
                    command['analyzed'] = False
                    
            
            measure_queue.publish(commands) # Send results for analysis






try:
    measure_queue
except NameError:
    # Only generate the measure_queue object if it
    # is not already defined.
    from CustomQueue import *
    measure_queue = Queue_measure()
