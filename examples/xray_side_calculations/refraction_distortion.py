#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Imports
########################################

import sys, os
SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

import glob
from SciAnalysis import tools
from SciAnalysis.Data import *
#from SciAnalysis.XSAnalysis.Data import *
#from SciAnalysis.XSAnalysis import Protocols



root_dir = './'
source_dir = os.path.join(root_dir, './')
output_dir = os.path.join(root_dir, './')






theta_incident = np.linspace(0, 1.0, num=5000)




# At 13.5 keV
# Ge: 0.186
#critical_angle = 0.186
# Si: 0.132
#critical_angle = 0.132
# SiO2: 0.141
#critical_angle = 0.141
# Polystyrene: 0.090    
#critical_angle = 0.090



critical_angle_film = 0.090
critical_angle_substrate = 0.132
critical_angle = critical_angle_substrate

lambda_A = 0.9184
k = 2.*np.pi/lambda_A
kpre = 2.*k
# Angles are full angles (2theta)
def q_to_angle(q):
    return np.degrees( 2.0*np.arcsin(q/kpre) )
def angle_to_q(angle):
    return kpre*np.sin(np.radians(angle/2))


if False:
    # Test nonlinear transformation
    lambda_A = 0.9184/100
    k = 2.*np.pi/lambda_A
    kpre = 2.*k
    # Angles are full angles (2theta)

    def q_to_angle(q):
        return np.degrees( 2.0*np.sqrt(np.arcsin(q/kpre)) )
    def angle_to_q(angle):
        return kpre*np.square(np.sin(np.radians(angle/2)))

class DataLines_current(DataLines):
    
    
    def _plot(self, save=None, show=False, plot_range=[None,None,None,None], plot_buffers=[0.2,0.05,0.2,0.05], error=False, error_band=False, xlog=False, ylog=False, xticks=None, yticks=None, dashes=None, **kwargs):
        
        plot_args = self.plot_args.copy()
        plot_args.update(kwargs)
        self.process_plot_args(**plot_args)
        
        self.fig = plt.figure( figsize=(10,10), facecolor='white' )
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        #fig_height = 1.0-top_buf-bottom_buf
        fig_height = fig_width
        self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )
        
        
        v_buf = 0.05
        fig_height_top = 1.0-top_buf-bottom_buf - v_buf - fig_height
        self.ax_top = self.fig.add_axes( [left_buf, bottom_buf+fig_height+v_buf, fig_width, fig_height_top] )
        
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
        
        
        # Axis scaling
        xi, xf, yi, yf = self.ax.axis()
        if plot_range[0] != None: xi = plot_range[0]
        if plot_range[1] != None: xf = plot_range[1]
        if plot_range[2] != None: yi = plot_range[2]
        if plot_range[3] != None: yf = plot_range[3]
        self.ax.axis( [xi, xf, yi, yf] )
        
        self._plot_extra(**plot_args)
        self._plot_extra_top(**plot_args)
        
        if save:
            if 'dpi' in plot_args:
                plt.savefig(save, dpi=plot_args['dpi'], transparent=True)
            else:
                plt.savefig(save)
        
        if show:
            self._plot_interact()
            plt.show()
            
        plt.close(self.fig.number)
            
            
    def _plot_main(self, error=False, error_band=False, dashes=None, **plot_args):
        
        for line in self.lines:
            
            
            
            plot_args_current = {}
            plot_args_current.update(self.plot_args)
            plot_args_current.update(plot_args)
            plot_args_current.update(line.plot_args)
            
            p_args = dict([(i, plot_args_current[i]) for i in self.plot_valid_keys if i in plot_args_current])
        

                
            if line.name=='critical_angle_exp':
                l, = self.ax.plot(line.x, line.y, label=line.name, markerfacecolor='#bfbf00', markeredgecolor='k', zorder=20, markeredgewidth=1, **p_args)
            elif line.name=='direct_beam_film':
                l, = self.ax.plot(line.x, line.y, label=line.name, **p_args)
                l.set_dashes([5,5])
            else:
                l, = self.ax.plot(line.x, line.y, label=line.name, **p_args)
                
            if dashes is not None:
                l.set_dashes(dashes)        
                            
    
    def _plot_extra(self, **plot_args):
        self.ax.set_aspect('auto')
        
        #two_theta_s_deg = q_to_angle(self.y_axis) # qz -> 2theta_s
        theta_deg = self.x_axis # theta_i
        
        # Horizon
        qzi = angle_to_q(theta_deg)
        self.ax.plot( self.x_axis, qzi, '-', linewidth=4.0, color='0.5', dashes=[15,15], zorder=-5)
        # Specular
        qzi = angle_to_q(theta_deg*2)
        self.ax.plot( self.x_axis, qzi, '-', linewidth=4.0, color='r', )
        # Yoneda film
        qzi = angle_to_q(theta_deg+critical_angle_film)
        self.ax.plot( self.x_axis, qzi, '-', linewidth=3.0, color='#bfbf00', zorder=-20, alpha=0.5)
        # Yoneda substrate
        qzi = angle_to_q(theta_deg+critical_angle_substrate)
        self.ax.plot( self.x_axis, qzi, '-', linewidth=3.0, color='#bfbf00', )
        
        

        xi, xf, yi, yf = self.ax.axis()
        self.ax.axis( [xi, xf, angle_to_q(xi), angle_to_q(xf)] )
        #self.ax.axis( [xi, xf, angle_to_q(xf)*-0.5, angle_to_q(xf)] )
        xi, xf, yi, yf = self.ax.axis()


        # Critical angle
        self.ax.axvline( critical_angle_film, color='0.75', linewidth=1.0, zorder=-20)
        self.ax.text(critical_angle_film, yf, r'$\theta_{cf}$', size=15, color='0.75', verticalalignment='bottom', horizontalalignment='center', zorder=-10)

        self.ax.axvline( critical_angle_substrate, color='0.75', linewidth=1.0, zorder=-20)
        self.ax.text(critical_angle_substrate, yf, r'$\theta_{cs}$', size=15, color='0.75', verticalalignment='bottom', horizontalalignment='center', zorder=-10)

        
        self.ax.text(q_to_angle(yf)*0.5, yf, r'$\mathrm{R}$', size=30, color='r', verticalalignment='bottom', horizontalalignment='center', zorder=-10)

        self.ax.text(q_to_angle(yf)-critical_angle_substrate, yf, r'$\mathrm{Y}$', size=30, color='#bfbf00', verticalalignment='bottom', horizontalalignment='center', zorder=-10)
        #self.ax.text(q_to_angle(yf)-critical_angle_film, yf, r'$\mathrm{Y}_{f}$', size=30, color='#bfbf00', verticalalignment='bottom', horizontalalignment='center', zorder=-20, alpha=0.5)

        self.ax.text(q_to_angle(yf), yf, r'$\mathrm{H}$', size=30, color='0.5', verticalalignment='bottom', horizontalalignment='center', zorder=-10)


        self.ax.text(xf, angle_to_q(0.02), r'$\mathrm{T}$', size=30, color='b', verticalalignment='bottom', horizontalalignment='right', zorder=-10)
        
        
        
        if True:
            self.ax2 = self.ax.twinx()

            yticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
            if False:
                # Assume linearity between axes
                #self.ax2.plot( self.x_axis, self.x_axis, 'b')
                self.ax2.axis( [xi, xf, q_to_angle(yi), q_to_angle(yf)] )
                self.ax2.set_yticks(yticks)
            else:
                # Use nonlinear mapping
                #self.ax2.plot( self.x_axis, angle_to_q(self.x_axis), 'b')
                
                yticks_transformed = []
                labels = []
                for tick in yticks:
                    qzpos = angle_to_q(tick)
                    yticks_transformed.append(qzpos)
                    
                    if tick<0:
                        tick = '−{}'.format(abs(tick))
                    elif tick==0:
                        tick = '{:.1f}'.format(abs(tick))
                    labels.append(tick)
                
                self.ax2.set_yticks(yticks_transformed)
                self.ax2.set_yticklabels(labels, size=25)
                
                self.ax2.set_ylabel(r'$2 \theta_s \, (^{\circ})$', size=30)
                
                self.ax2.axis( [xi, xf, yi, yf] )
                
                
    def _plot_extra_top(self, **plot_args):
        
        line = self.extra_line
        self.ax_top.plot(line.x, line.y, '-', linewidth=4.0, color='0.5', )
        
        self.ax_top.set_xlabel('')
        self.ax_top.set_xticks(self.ax.get_xticks())
        self.ax_top.set_xticklabels([])
        self.ax_top.set_ylabel('$\Delta_{q_z}$', size=20)
        for tick in self.ax_top.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)
            
                
        

lines = DataLines_current()
lines.x_axis = theta_incident

if False:
    # Horizon
    
    line = DataLine(x=theta_incident, y=angle_to_q(horizon))
    lines.add_line(line)

    line.plot_args['color'] = '0.5'
    line.plot_args['linestyle'] = '-'
    line.plot_args['linewidth'] = 4.0
    line.plot_args['marker'] = None
    line.plot_args['dashes'] = [5,5]

if False:
    # Specular reflection
    specular = 2.0*theta_incident
    line = DataLine(x=theta_incident, y=angle_to_q(specular))
    lines.add_line(line)

    line.plot_args['color'] = 'purple'
    line.plot_args['linestyle'] = '-'
    line.plot_args['linewidth'] = 4.0
    line.plot_args['marker'] = None


if False:
    incident_angle_exp = npzfile['incident_angle_experiment']
    critical_angle_exp = npzfile['Yoneda_Ge_experiment']
    #print(incident_angle_exp)


    line = DataLine(x=incident_angle_exp, y=critical_angle_exp, name='critical_angle_expJL')
    lines.add_line(line)

    line.plot_args['color'] = 'purple'
    line.plot_args['linestyle'] = 'o'
    line.plot_args['linewidth'] = 4.0
    line.plot_args['markersize'] = 12.0





    

horizon = 1.0*theta_incident
specular = 2.0*theta_incident

if True:
    # T
    alpha_i_rad = np.arccos( np.cos(np.radians(theta_incident))/np.cos(np.radians(critical_angle_substrate)) )
    alpha_i_deg = np.degrees(alpha_i_rad)
    two_theta_s_deg = np.where( theta_incident>critical_angle_substrate, theta_incident - alpha_i_deg, horizon )
    T = angle_to_q(two_theta_s_deg)


    line = DataLine(x=theta_incident, y=T, name='direct_beam_substrate')
    lines.add_line(line)

    line.plot_args['color'] = 'b'
    line.plot_args['linestyle'] = '-'
    line.plot_args['linewidth'] = 4.0
    line.plot_args['marker'] = None
    #line.plot_args['markersize'] = 12.0



if True:
    # T
    alpha_i_rad = np.arccos( np.cos(np.radians(theta_incident))/np.cos(np.radians(critical_angle_film)) )
    alpha_i_deg = np.degrees(alpha_i_rad)
    two_theta_s_deg = np.where( theta_incident>critical_angle_film, theta_incident - alpha_i_deg, horizon )
    T = angle_to_q(two_theta_s_deg)


    line = DataLine(x=theta_incident, y=T, name='direct_beam_film')
    lines.add_line(line)

    line.plot_args['color'] = 'b'
    line.plot_args['linestyle'] = '-'
    line.plot_args['linewidth'] = 4.0
    line.plot_args['marker'] = None
    #line.plot_args['markersize'] = 12.0





# Show the detector position (qz) of candidate reciprocal-space (Qz) peaks

#for Qz_consider in np.arange(0.0, 0.20, 0.005):
#for Qz_consider in np.arange(-0.2, 0.20, 0.01):
#for Qz_consider in [0.01, 0.03, 0.06, 0.09]:
for Qz_consider in [0.05]:


    two_alpha_s_deg = q_to_angle(Qz_consider)
    #two_alpha_s_deg_vector = np.ones(len(theta_incident))*two_alpha_s_deg
    
    
    qz_measure_T = np.zeros(len(theta_incident))
    qz_measure_R = np.zeros(len(theta_incident))


    
    # True position
    if True:
        Qz_true = np.ones(len(theta_incident))*Qz_consider
        line = DataLine(x=theta_incident, y=Qz_true, name='Qz_true')
        lines.add_line(line)

        line.plot_args['color'] = 'blue'
        line.plot_args['linestyle'] = '-'
        line.plot_args['linewidth'] = 1.0
        line.plot_args['marker'] = None
        line.plot_args['alpha'] = 0.2
                
    
    # T channel refracted peak
    if True:
        
        # Above horizon component
        if True:
            
            alpha_i_rad = np.arccos( np.cos(np.radians(theta_incident))/np.cos(np.radians(critical_angle_film)) )
            alpha_i_deg = np.where( theta_incident>critical_angle_film, np.degrees(alpha_i_rad), 0.0 )
            
            # Above-horizon
            alpha_f_deg = two_alpha_s_deg - alpha_i_deg
            alpha_f_rad = np.radians(alpha_f_deg)
            
            theta_f_rad = np.arccos( np.cos(np.radians(critical_angle_film))*np.cos(alpha_f_rad) )
            theta_f_deg = np.where( alpha_i_deg<two_alpha_s_deg, np.degrees(theta_f_rad), np.nan ) # Only valid above horizon
            
            
            two_theta_s_deg = theta_incident + theta_f_deg
            two_theta_s_deg = np.where( theta_incident>critical_angle_film, two_theta_s_deg, np.nan)
            qz = angle_to_q(two_theta_s_deg)
            
            qz_measure_T += np.where(np.isnan(qz), 0, qz)
            
            line = DataLine(x=theta_incident, y=qz, name='Tpeak')
            lines.add_line(line)

            line.plot_args['color'] = 'blue'
            line.plot_args['linestyle'] = '-'
            line.plot_args['linewidth'] = 3.0
            line.plot_args['marker'] = None
            line.plot_args['alpha'] = 0.8
            
            # Below film critical angle
            # (incident beam traveling along horizon)
            if True:
                two_theta_s_deg = theta_incident + theta_f_deg
                two_theta_s_deg = np.where( theta_incident<critical_angle_film, two_theta_s_deg, np.nan)
                qz = angle_to_q(two_theta_s_deg)
                
                qz_measure_T += np.where(np.isnan(qz), 0, qz)
                qz_measure_R += np.where(np.isnan(qz), 0, qz)
                
                line = DataLine(x=theta_incident, y=qz, name='Tpeak')
                lines.add_line(line)

                line.plot_args['color'] = 'blue'
                line.plot_args['linestyle'] = '-'
                line.plot_args['linewidth'] = 3.0
                line.plot_args['marker'] = None    
                line.plot_args['alpha'] = 0.2
                        
        
        # Below-horizon (film refraction only)
        if False:
            
            two_theta_s_deg = np.where( alpha_i_deg>two_alpha_s_deg,  theta_incident - alpha_i_deg + two_alpha_s_deg , np.nan ) # Only valid below horizon
            qz = angle_to_q(two_theta_s_deg)
            
            line = DataLine(x=theta_incident, y=qz, name='Tpeak')
            lines.add_line(line)

            line.plot_args['color'] = 'blue'
            line.plot_args['linestyle'] = '-'
            line.plot_args['linewidth'] = 3.0
            line.plot_args['marker'] = None    
            line.plot_args['alpha'] = 0.2

        # Below-horizon (film and substrate refractions)
        if True:
            
            incident_prime_deg = alpha_i_deg - two_alpha_s_deg
            incident_prime_rad = np.radians(incident_prime_deg)
            # Refraction at film-substrate interface
            n_ratio = np.cos(np.radians(critical_angle_substrate))/np.cos(np.radians(critical_angle_film))
            output_prime_rad = np.arccos( np.cos(incident_prime_rad)/n_ratio )
            output_prime_deg = np.degrees(output_prime_rad)
            
            two_theta_s_deg = theta_incident - output_prime_deg
            two_theta_s_deg = np.where( alpha_i_deg>two_alpha_s_deg, two_theta_s_deg, np.nan ) # Only valid below horizon
            two_theta_s_deg = np.where( theta_incident>critical_angle_film, two_theta_s_deg, np.nan ) # Don't show when incident beam below critical angle
            qz = angle_to_q(two_theta_s_deg)
            
            qz_measure_T += np.where(np.isnan(qz), 0, qz)
            
            line = DataLine(x=theta_incident, y=qz, name='Tpeak')
            lines.add_line(line)

            line.plot_args['color'] = 'blue'
            line.plot_args['linestyle'] = '-'
            line.plot_args['linewidth'] = 3.0
            line.plot_args['marker'] = None    
            line.plot_args['alpha'] = 0.8
            
            # Scattering refracted along horizon
            if True:
                two_theta_s_deg = theta_incident - output_prime_deg
                two_theta_s_deg = np.where( np.isnan(two_theta_s_deg) , horizon, np.nan)
                two_theta_s_deg = np.where( alpha_i_deg>two_alpha_s_deg, two_theta_s_deg, np.nan ) # Only valid below horizon
                qz = angle_to_q(two_theta_s_deg)
                
                qz_measure_T += np.where(np.isnan(qz), 0, qz)

                line = DataLine(x=theta_incident, y=qz, name='Tpeak')
                lines.add_line(line)

                line.plot_args['color'] = 'blue'
                line.plot_args['linestyle'] = '-'
                line.plot_args['linewidth'] = 3.0
                line.plot_args['marker'] = None    
                line.plot_args['alpha'] = 0.2
                
                
            # Below film critical angle
            # (incident beam traveling along horizon)
            if True and two_alpha_s_deg<0:
                
                angle_vals = theta_incident - output_prime_deg
                angle_vals = np.where( theta_incident<critical_angle_film, angle_vals, np.nan )
                qz = angle_to_q(angle_vals)
                
                line = DataLine(x=theta_incident, y=qz, name='Tpeak')
                lines.add_line(line)

                line.plot_args['color'] = 'blue'
                line.plot_args['linestyle'] = '-'
                line.plot_args['linewidth'] = 3.0
                line.plot_args['marker'] = None    
                line.plot_args['alpha'] = 0.2    
                


            
            
    # R channel refracted peak
    if True:
        
        if True:
            alpha_i_rad = np.arccos( np.cos(np.radians(theta_incident))/np.cos(np.radians(critical_angle_film)) )
            alpha_i_deg = np.where( theta_incident>critical_angle_film, np.degrees(alpha_i_rad), 0.0 )
            
            alpha_f_rad = np.radians(two_alpha_s_deg) + alpha_i_rad
            theta_f_rad = np.arccos( np.cos(np.radians(critical_angle_film))*np.cos(alpha_f_rad) )
            theta_f_deg = np.degrees(theta_f_rad)

            two_theta_s_deg = theta_incident + theta_f_deg

            qz = angle_to_q(two_theta_s_deg)
            
            qz_measure_R += np.where(np.isnan(qz), 0, qz)
            
            line = DataLine(x=theta_incident, y=qz, name='Rpeak')
            lines.add_line(line)
        
            line.plot_args['color'] = 'red'
            line.plot_args['linestyle'] = '-'
            line.plot_args['linewidth'] = 3.0
            line.plot_args['marker'] = None
            line.plot_args['alpha'] = 0.8        
            
            
            
        # Scattering goes below horizon
        if True and two_alpha_s_deg<0:
            # For this to occur:
            #  two_alpha_s_deg < 0
            #  |two_alpha_s_deg| > |alpha_i_deg|
            
            incident_prime_deg = np.where( theta_incident>critical_angle_film, np.abs(two_alpha_s_deg) - alpha_i_deg, np.nan)
            incident_prime_rad = np.radians(incident_prime_deg)
            # Refraction at film-substrate interface
            n_ratio = np.cos(np.radians(critical_angle_substrate))/np.cos(np.radians(critical_angle_film))
            output_prime_rad = np.arccos( np.cos(incident_prime_rad)/n_ratio )
            output_prime_deg = np.degrees(output_prime_rad)
            
            two_theta_s_deg = np.where( np.abs(two_alpha_s_deg)>np.abs(alpha_i_deg), theta_incident - output_prime_deg, np.nan )
            
            qz = angle_to_q(two_theta_s_deg)
            
            line = DataLine(x=theta_incident, y=qz, name='Rpeak')
            lines.add_line(line)
        
            line.plot_args['color'] = 'red'
            line.plot_args['linestyle'] = '-'
            line.plot_args['linewidth'] = 3.0
            line.plot_args['marker'] = None
            line.plot_args['alpha'] = 0.8        
            
            
            # Scattering traveling along horizon
            # (at film-substrate interface, since beam is hitting substrate interface below its critical angle)
            if True:
                
                internal_critical_angle_rad = np.arccos(n_ratio)
                
                angle_vals = np.where( incident_prime_rad<internal_critical_angle_rad, horizon, np.nan )
                angle_vals = np.where( incident_prime_rad>0, angle_vals, np.nan )
                qz = angle_to_q(angle_vals)
                
                line = DataLine(x=theta_incident, y=qz, name='Rpeak')
                lines.add_line(line)
                
                line.plot_args['color'] = 'red'
                line.plot_args['linestyle'] = '-'
                line.plot_args['linewidth'] = 3.0
                line.plot_args['marker'] = None    
                line.plot_args['alpha'] = 0.2  
        
    
    
    
    # Amount of refraction distortion
    if True:
        
        Qz_true = np.ones(len(theta_incident))*Qz_consider
        
        Delta_qz_T = qz_measure_T - Qz_true

        Delta_qz_R = qz_measure_R - Qz_true
        
        delta_Delta = Delta_qz_R - Delta_qz_T
        
        delta_Delta_m = delta_Delta - angle_to_q(2*theta_incident)
        
        line = DataLine(x=theta_incident, y=delta_Delta, name='Delta')
        
        lines.extra_line = line

        if False:
            lines.add_line(line)
            line.plot_args['color'] = 'cyan'
            line.plot_args['linestyle'] = '-'
            line.plot_args['linewidth'] = 3.0
            line.plot_args['marker'] = None    
            line.plot_args['alpha'] = 0.8   







    #Delta_qz = qz - kpre*np.sin( 0.5*np.arccos(np.cos(theta_f_rad)/np.cos(np.radians(critical_angle_film))) + 0.5*np.arccos(np.cos(np.radians(theta_incident))/np.cos(np.radians(critical_angle_film))) )
    #qz = Qz_consider + Delta_qz

    #line = DataLine(x=theta_incident, y=qz, name='Delta')
    #lines.add_line(line)

    #line.plot_args['color'] = 'purple'
    #line.plot_args['linestyle'] = '-'
    #line.plot_args['linewidth'] = 2.0
    #line.plot_args['marker'] = None
    
    
    
    #test_f_rad = np.arccos( np.cos(np.radians(critical_angle_film))*np.cos(two_alpha_s_deg_vector) )
    #test_f_deg = np.degrees(test_f_rad)


    #line = DataLine(x=theta_incident, y=angle_to_q(alpha_f_deg), name='show')
    #lines.add_line(line)

    #line.plot_args['color'] = 'cyan'
    #line.plot_args['linestyle'] = '-'
    #line.plot_args['linewidth'] = 2.0
    #line.plot_args['marker'] = None


    

lines.plot_args = { 'color' : 'k',
                'marker' : 'o',
                'linewidth' : 3.0,
                'legend_location': 'NE'
                }  

lines.plot_args['rcParams'] = { 
                'axes.labelsize': 45,
                'xtick.labelsize': 35,
                'ytick.labelsize': 35,    
                #'legend.borderpad' : 0 ,
                'legend.fontsize' : 30 ,
                #'legend.numpoints' : 1 ,
                #'legend.handlelength' : 1.0 ,
                'legend.labelspacing' : 0.25 ,
                'legend.handletextpad' : 0.5 ,
                #'legend.columnspacing' : 0.0 ,
                'xtick.major.pad': 14,
                'ytick.major.pad': 14,
                
                }
    

lines.x_label = 'angle'
lines.x_rlabel = r'$\theta_i \, (^{\circ})$'
lines.y_label = 'qz'
lines.y_rlabel = '$q_z \, (\mathrm{\AA^{-1}})$'
outfile = os.path.join(output_dir, 'fig-refraction_distortion.png')
lines.plot(save=outfile, plot_range=[0, 1.0, None, None], plot_buffers=[0.27, 0.16, 0.16, 0.05], _xticks=[0, 0.04, 0.08, 0.12], dpi=200)
    




