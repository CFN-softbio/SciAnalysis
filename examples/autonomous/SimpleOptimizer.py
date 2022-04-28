#!/usr/bin/python3
import numpy as np

class SimpleOptimizer():
    
    def __init__(self, parameter_ranges, initial_guess=None, evaluation_function=None, name='SimpOpt', verbosity=3, **kwargs):
        
        self.kwargs = kwargs
        self.parameter_ranges = np.asarray(parameter_ranges)
        self.num_parameters = len(self.parameter_ranges)
        
        if initial_guess is None:
            initial_guess = [ 0.5*(lower+upper) for lower, upper in parameter_ranges ]
        self.initial_guess = np.asarray(initial_guess)
        if evaluation_function is None:
            evaluation_function = self.evaluate
        self.evaluation_function = evaluation_function
        
        self.name = name
        self.verbosity = verbosity
        
        self.pos_best, self.err_best = None, None
        
    # Helpers
    ########################################
    def msg(self, txt, threshold=3, indent=0, indent_txt='  ', verbosity=None, empty_lines=0, **kwargs):
        '''Outputs a status line indicating the current state of execution.'''
        if verbosity is None:
            verbosity = self.verbosity
        if verbosity>=threshold:
            indent = np.clip(indent, 0, 10)
            indent = indent_txt*indent
            for i in range(empty_lines):
                print('')
            print('{}> {}{}'.format(self.name, indent, txt))        
        

    def best(self):
        if self.pos_best is None:
            self.pos_best = self.initial_guess
        if self.err_best is None:
            self.err_best = self.evaluation_function(self.pos_best)
        return self.pos_best, self.err_best


    # Optimization algorithms
    ########################################
        
    def evaluate(self, position):
        '''Test evaluation function, which simply wants the position to be zero.'''
        err = 0
        for p in position:
            err += np.abs(p)
            
        return err


    def optimize_not(self):
        '''Shell that can be used as a starting point.'''
        
        pos_best, err_best = self.best()

        self.pos_best, self.err_best = pos_best, err_best
        return pos_best, err_best
    
    
    def optimize_random(self, iterations=10, **kwargs):
        
        pos_best, err_best = self.best()
        
        for i in range(iterations):
            pos = np.asarray([ np.random.uniform(lower, upper) for lower, upper in self.parameter_ranges ])
            err = self.evaluation_function(pos)
            if err<err_best:
                pos_best = pos
                err_best = err
            self.msg('random {}\t{:,.1f}\t{}'.format(i, err_best, pos_best), 4, 1)
                
        self.pos_best, self.err_best = pos_best, err_best
        return pos_best, err_best
        
        
    def optimize_local_step(self, iterations=20, stride_rel_initial=0.04, stride_rel_update=0.2, **kwargs):
        '''Naive local stepping (with progressively smaller step size) in the 
        downward direction. This only works for very well-behaved error functions
        (or else walks down in the local well).'''

        pos_best, err_best = self.best()
        
        spans = np.asarray([ abs(upper-lower) for lower, upper in self.parameter_ranges ])
        step_sizes = spans*stride_rel_initial
        polarities = np.ones(self.num_parameters)
        
        pos = np.copy(pos_best)
        i = 0
        while(i<iterations):
            
            for ip in range(self.num_parameters):
                fails = 0
                count_steps = 0
                while(fails<2):
                    if count_steps>5:
                        # If we keep walking in this direction
                        # (it's productive since it keeps lowering error),
                        # we can start taking bigger steps
                        step_sizes[ip] *= 1.25/stride_rel_update
                        count_steps = 0
                    
                    # Take a candidate step
                    pos_new = np.copy(pos)
                    pos_new[ip] += polarities[ip]*step_sizes[ip]
                    
                    if self.verbosity>=10:
                        txt = '    parameter{}: step {:+.2g} to [{:,.4g}]'.format(ip, polarities[ip]*step_sizes[ip], pos_new[ip])
                    if pos_new[ip]>self.parameter_ranges[ip][0] and pos_new[ip]<self.parameter_ranges[ip][1]:
                        # This point is within range
                        err = self.evaluation_function(pos_new)
                        count_steps += 1
                    
                        if err<err_best:
                            if self.verbosity>=10:
                                print('{}; is better ({:,.6g})'.format(txt, err))
                            # This point is indeed better
                            pos = pos_new
                            pos_best = pos
                            err_best = err
                        else:
                            if self.verbosity>=10:
                                print('{}; is worse ({:,.6g})'.format(txt, err))
                            # Ignore this movement; go in opposite direction
                            polarities[ip] *= -1
                            step_sizes[ip] *= stride_rel_update
                            fails += 1
                    else:
                        if self.verbosity>=10:
                            print('{}; hits boundary'.format(txt))
                        # This point is out-of-range
                        polarities[ip] *= -1
                        fails += 1
                        
            self.msg('lw {}\t{:,.1f}\t{}'.format(i, err_best, pos_best), 4, 1)
            i += 1
        
        self.pos_best, self.err_best = pos_best, err_best
        return pos_best, err_best
        
        
    def optimize_siman(self, T_initial=90, T_final=0.1, alpha=1.0, stride_rel_initial=0.1, stride_rel_update=0.94, e_normalization=None, **kwargs):
        '''Simple implementation of simulated annealing.'''
        
        pos_best, err_best = self.best()
        if e_normalization is None:
            e_normalization = 1.0
        elif e_normalization=='auto':
            e_normalization = err_best/100
        
        spans = np.asarray([ abs(upper-lower) for lower, upper in self.parameter_ranges ])
        step_sizes = spans*stride_rel_initial
        step_rel = stride_rel_initial
        
        T_current = T_initial
        pos_c = pos_best
        err_c = err_best
        
        i = 0
        while T_current>T_final:
            
            pos_new = np.copy(pos_c)
            valid = True
            for ip in range(self.num_parameters):
                pos_new[ip] += np.random.uniform(-step_sizes[ip], +step_sizes[ip])
                if pos_new[ip]<self.parameter_ranges[ip][0] or pos_new[ip]>self.parameter_ranges[ip][1]:
                    valid = False
                    
            if valid: # Ignore points that put us out of bounds
                err_new = self.evaluation_function(pos_new)
                
                if err_new<err_c:
                    # Accept improved position
                    pos_c = pos_new
                    err_c = err_new
                else:
                    # Accept 'worse' position with a probability
                    if np.random.uniform(0,1) < np.exp( -(err_new-err_c)/(T_current*e_normalization) ):
                        pos_c = pos_new
                        err_c = err_new
                        
                if err_new<err_best:
                    pos_best = pos_new
                    err_best = err_new
                
                    
            T_current -= alpha
            step_sizes *= stride_rel_update
            step_rel *= stride_rel_update
            self.msg('siman {} (T={:.1f}, step={:.2g})\t{:,.1f}\t{}'.format(i, T_current, step_rel, err_c, pos_c), 4, 1)
            i += 1
                    
                
        self.pos_best, self.err_best = pos_best, err_best
        return pos_best, err_best
        
        
        
        
        
        
        
        
        
        
        
