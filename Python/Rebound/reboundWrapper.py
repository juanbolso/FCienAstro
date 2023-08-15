import rebound
import numpy as np
import matplotlib.pyplot as plt
from math import *
from ipywidgets import IntProgress
# ***** NUMERICAL INTEGRATIONS *****


# Conversion factors:
G2R = pi/180
R2G = 180/pi


class Integration2coplanarPlanets:

    def __init__(self, m0, m1, m2, a1i, e1i, w1i, M1i, a2i, e2i, w2i, M2i, tTot, dT, dOut, RHtol = 2*np.sqrt(3), MMRvariables = True, K1=1, K2=2):
            
        # Create simulation
        self.sim = rebound.Simulation()

        # Integration algorithm:
        # self.sim.integrator = "whfast"
        self.sim.integrator = "mercurius" 
        self.sim.ri_mercurius.hillfac = RHtol # Hill radius for mercurius

        # Reference system:
        # self.sim.ri_whfast.coordinates = 'jacobi' #default
        # self.sim.ri_whfast.coordinates = 'democraticheliocentric' # ~ Poincare
        # self.sim.ri_whfast.coordinates = 'whds'
        self.sim.ri_mercurius.coordinates = 'democraticheliocentric' # ~ Poincare

        # Units and time parameters:
        self.sim.units = ('yrs', 'AU', 'Msun')   
        #self.sim.t = tTot          # total int. time
        self.sim.dt = dT            # paso de la int.
        self.Nout  = int(tTot/dOut) # Total int. points = total int. time / output step
        # Improve performance & accuaracy:
        # self.sim.ri_whfast.safe_mode = 0
        # self.sim.ri_whfast.corrector = 11 # Solo en Jacobi
        KG2 = sqrt(self.sim.G)

        # DOBJECTS DEFINITION:
        self.sim.add(m=m0)                                                                  # Star        
        self.sim.add(m=m1/1047.57, a=a1i, e=e1i, inc=0*G2R, pomega=w1i*G2R, Omega=0*G2R, M=M1i*G2R) # Inner planet
        self.sim.add(m=m2/1047.57, a=a2i, e=e2i, inc=0*G2R, pomega=w2i*G2R, Omega=0*G2R, M=M2i*G2R) # Outer planet
        self.particles = self.sim.particles
        #self.sim.status()

        # Variable initiallization:
        self.times = np.linspace(0, tTot, self.Nout)
        self.t1 = self.times
        self.t2 = self.times
        self.a1 = np.zeros(self.Nout)
        self.e1 = np.zeros(self.Nout)
        self.w1 = np.zeros(self.Nout)
        self.l1 = np.zeros(self.Nout)
        self.a2 = np.zeros(self.Nout)
        self.e2 = np.zeros(self.Nout)
        self.w2 = np.zeros(self.Nout)
        self.l2 = np.zeros(self.Nout)
        self.MMRvariables = MMRvariables
        if self.MMRvariables:
            # MMR quantitites:
            self.K1 = K1
            self.K2 = K2
            self.deltaW = np.zeros(self.Nout)
            self.theta = np.zeros(self.Nout)
            self.s1 = np.zeros(self.Nout)
            self.s2 = np.zeros(self.Nout)
            self.ds = np.zeros(self.Nout)

        # Initiallize some possible used variables
        self.angle_yticks = None
        self.figure1 = None
        self.figure2 = None
        self.figure3 = None


    def startIntegration(self):
        # Progress bar indicator:
        print('Progress:')
        barra = IntProgress(min=0, max=self.Nout) # instantiate the bar
        display(barra) # display the bar

        # *** INTEGRATION STARTS ***
        # Integrate in steps (Ttot/Nout) y save the elements:
        for i,time in enumerate(self.times):
            self.sim.integrate(time, exact_finish_time=0)
            self.a1[i] = self.particles[1].a
            self.e1[i] = self.particles[1].e
            self.w1[i] = (self.particles[1].pomega)*R2G % 360
            self.l1[i] = (self.particles[1].l)*R2G % 360
            self.a2[i] = self.particles[2].a
            self.e2[i] = self.particles[2].e
            self.w2[i] = (self.particles[2].pomega)*R2G % 360
            self.l2[i] = (self.particles[2].l)*R2G % 360
            barra.value += 1 # signal to increment the progress bar
        if self.MMRvariables:
            # MMR quantitites:
            self.deltaW = (self.w1-self.w2) % 360
            self.theta = (self.K1*self.l1-self.K2*self.l2) % 360
            self.s1 = (self.theta + (self.K2-self.K1)*self.w1) % 360
            self.s2 = (self.theta + (self.K2-self.K1)*self.w2) % 360
            self.ds = (self.theta + (self.K2-self.K1)*self.deltaW) % 360
        #sim.status()

        # Just in case:
        del self.sim.particles
        self.sim = None # free the memory
        # *** INTEGRATION ENDS ***

        # Calculate mean a1 to correct a1i (structure law)
        #print(a1.mean())


    def proccessAngles(self, positiveAngles = True):
        # Angles shifting:
        self.w1 = self.w1 % 360
        self.w2 = self.w2 % 360
        self.l1 = self.l1 % 360
        self.l2 = self.l2 % 360
        if not positiveAngles:
            self.w1 = np.where(self.w1>180, self.w1 - 360, self.w1)
            self.l1 = np.where(self.l1>180, self.l1 - 360, self.l1)
            self.s1 = np.where(self.s1>180, self.s1 - 360, self.s1)
            self.w2 = np.where(self.w2>180, self.w2 - 360, self.w2)
            self.l2 = np.where(self.l2>180, self.l2 - 360, self.l2)
            self.s2 = np.where(self.s2>180, self.s2 - 360, self.s2)
            self.deltaW = np.where(self.deltaW>180, self.deltaW - 360, self.deltaW)
            self.theta = np.where(self.theta>180, self.theta - 360, self.theta)
            self.ds = np.where(self.ds>180, self.ds - 360, self.ds)
            self.angle_yticks = np.arange(-180,270,90) 
        else:
            self.angle_yticks = np.arange(0,450,90)


    def plotIntegration(self, size = (16, 9)):
        # *** START PLOTS ***

        # PLANET 1 PLOTS
        self.figure1, axs = plt.subplots(4, 1, figsize = (size[0], size[1]), sharex = 'col', sharey = False)
        # Gráficos: a(t), e(t) y ángulos(t).
        axs[0].plot(self.t1, self.a1, 's', marker = '.', ms = 1, c='k')
        axs[1].plot(self.t1, self.e1, 's', marker = '.', ms = 1, c='k')
        axs[2].plot(self.t1, self.w1, 's', marker = '.', ms = 1, c='k')
        axs[3].plot(self.t1, self.s1, 's', marker = '.', ms = 1, c='k')
        # Labels y ticks:
        axs[0].set(ylabel = '$a_1$ (au)', ylim = [0.9995*min(self.a1), 1.0005*max(self.a1)])
        axs[1].set(ylabel = '$e_1$', ylim = [0, 1], yticks = [0, 0.5 ,1])
        axs[2].set(ylabel = r'$\varpi_1$ (°)', yticks = self.angle_yticks)
        axs[3].set(ylabel = r'$\sigma_1$ (°)', xlabel = 't (years)', yticks = self.angle_yticks)
        plt.show()

        # PLANET 2 PLOTS
        self.figure2, axs = plt.subplots(4, 1, figsize = (size[0], size[1]), sharex = 'col', sharey = False)
        # Gráficos: a(t), e(t) y ángulos(t).
        axs[0].plot(self.t2, self.a2, 's', marker = '.', ms = 1, c='k')
        axs[1].plot(self.t2, self.e2, 's', marker = '.', ms = 1, c='k')
        axs[2].plot(self.t2, self.w2, 's', marker = '.', ms = 1, c='k')
        axs[3].plot(self.t2, self.s2, 's', marker = '.', ms = 1, c='k')
        # Labels y ticks:
        axs[0].set(ylabel = '$a_2$ (au)', ylim = [0.9995*min(self.a2), 1.0005*max(self.a2)])
        axs[1].set(ylabel = '$e_2$', ylim = [0, 1], yticks = [0, 0.5 ,1])
        axs[2].set(ylabel = r'$\varpi_2$ (°)', yticks = self.angle_yticks)
        axs[3].set(ylabel = r'$\sigma_2$ (°)', xlabel = 't (years)', yticks = self.angle_yticks)
        plt.show()

        # OTHER PLOTS
        self.figure3, axs = plt.subplots(3, 1, figsize = (size[0], size[1]), sharex = 'col', sharey = False)
        # Gráficos: a(t), e(t) y ángulos(t).
        axs[0].plot(self.t1, self.deltaW, 's', marker = '.', ms = 1, c='k')
        axs[1].plot(self.t1, self.theta, 's', marker = '.', ms = 1, c='k')
        axs[2].plot(self.t1, self.ds, 's', marker = '.', ms = 1, c='k')
        # Labels y ticks:
        axs[0].set(ylabel = r'$\Delta\varpi$ (°)', yticks = self.angle_yticks)
        axs[1].set(ylabel = r'$\theta$ (°)', yticks = self.angle_yticks)
        axs[2].set(ylabel = r'$\sigma^*$ (°)', xlabel = 't (years)', yticks = self.angle_yticks);
        plt.show()


    def savePlots(self, dir, filePlot1, filePlot2, filePlot3, format = "png"):
        if not os.path.exists(dir):
            os.makedirs(dir)

        self.figure1.savefig(dir + filePlot1, bbox_inches = 'tight', dpi = int(cal), format = format);  
        print('Plot 1 saved in:', dir+ filePlot1 + '.' + format)

        self.figure2.savefig(dir + filePlot2, bbox_inches = 'tight', dpi = int(cal), format = format);  
        print('Plot 2 saved in:', dir + filePlot2 + '.' + format)

        self.figure3.savefig(dir + filePlot3, bbox_inches = 'tight', dpi = int(cal), format = format);  
        print('Plot 3 saved in:', dir + filePlot3 + '.' + format)


    def saveData(self, dir, fileData1, fileData2):

        # data format and headers
        formato_data = "%6.2f %8.6f %8.6f %5.2f %5.2f"
        nom_cols1 = "    t    a1       e1      w1  lambda1"
        nom_cols2 = "    t    a2       e2      w2  lambda2"

        # SAVE NUM INT DATA
        data1 = np.c_[self.t1, self.a1, self.e1, self.w1, self.l1]
        data2 = np.c_[self.t2, self.a2, self.e2, self.w2, self.l2]
        np.savetxt(dir + fileData1, data1, fmt=formato_data, header = nom_cols1, comments='')
        np.savetxt(dir + fileData2, data2, fmt=formato_data, header = nom_cols2, comments='')
        print('Data saved in:', dir)

