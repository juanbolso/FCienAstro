import pandas as pd
import numpy as np
import scipy as sc
from math import *


# *** HELPER FUNCTIONS ***
G2R = pi/180 

# FUNCIÓN PARA OBTENER EL ANGULO MAS CERCANO DE UN ARRAY DE ANGULOS:
def find_nearest_ang(ang_arr, ang):
    ind = 0
    angulo = ang_arr[ind]
    X = sin(ang*G2R)
    Y = cos(ang*G2R)
    deltamin = 2    
    for i in range(len(ang_arr)):    
        angi = ang_arr[i]    
        x = sin(angi*G2R)
        y = cos(angi*G2R)    
        delta = (x-X)**2 + (y-Y)**2   
        if (delta < deltamin):
            deltamin = delta
            ind = i
            angulo = angi   
    return ind
    #return angulo, ind

def avg_ang(a, b): # https://en.wikipedia.org/wiki/Circular_mean
    ax = sin(a*G2R)
    bx = sin(b*G2R)
    ay = cos(a*G2R)
    by = cos(b*G2R)
    avg = (atan2((ax+bx)/2,(ay+by)/2))/G2R
    return avg

def interpolateArray(X, Y):
    num_interpolated_values = len(X) - len(Y) 
    step_size = num_interpolated_values // (len(Y) - 1)
    interpolated_Y = []

    for i in range(len(Y) - 1):
        start_value = Y[i]
        end_value = Y[i + 1]
        interpolated_segment = np.linspace(start_value, end_value, step_size + 2)[1:-1]
        interpolated_Y.extend(interpolated_segment)

    interpolated_Y.append(Y[-1])
    return interpolated_Y



class Hsurface:

    def __init__(self, dirH: str, fileH1: str, fileH2: str, K1: int, K2: int, positiveAngles = True, mapMindis = False, removeOutliers = True, RHtol = 2*np.sqrt(3)):
        # ATRIBUTES
        self.K1 = K1
        self.K2 = K2
        self.mapMindis = mapMindis
        #self.mindis = []
        # LOAD & MERGE H SURFACES DATA:
        self.dataHsup = pd.merge(self.loadHsurface(dirH + fileH1), self.loadHsurface(dirH + fileH2), how='outer')
        # IMPROVE H SURFACE DATA:
        self.improveHdata(removeOutliers, RHtol, 0.99)
        # Obtain variables from merged dataframe:
        self.e1 = self.dataHsup['e1'].values
        self.e2 = self.dataHsup['e2'].values
        self.deltaW = self.dataHsup['deltavarpi'].values % 360
        self.sigma = self.dataHsup['sigma1'].values % 360
        self.H = -self.dataHsup['R'].values # Here it could be added the H_0 (Keplerian) term in the future
        self.theta = (self.sigma - (self.K2-self.K1)*self.deltaW) % 360
        # MAP MINDIS (OR DELTA) SURFACE?
        if self.mapMindis: self.mindis = self.dataHsup['mindis'].values
        # ANGLES DOMAIN ADEQUATION:
        self.positiveAngles = positiveAngles
        self.anglesDomain(self.positiveAngles)
        self.reOrderColumns()

    def loadHsurface(self, path):
        return pd.read_csv(path, delimiter=r"\s+")

    def improveHdata(self,removeOutliers, RHtol, eLimit):
      # To avoid potentially e=1 noise, we filter those value
        self.dataHsup = self.dataHsup[self.dataHsup['e1'].values < eLimit] 
        self.dataHsup = self.dataHsup[self.dataHsup['e2'].values < eLimit] 
        # Remove outliers? (https://www.kite.com/python/answers/how-to-remove-outliers-from-a-pandas-dataframe-in-python) 
        H = self.dataHsup['R'].values
        if removeOutliers: self.dataHsup = self.dataHsup[(np.abs(sc.stats.zscore(H)) < 3)] 
        # Filter points in encounter condition according RHtol
        self.dataHsup = self.dataHsup[self.dataHsup['mindis'].values >= RHtol]

    def anglesDomain(self, positiveAngles):
        if not positiveAngles: 
            self.deltaW = np.where(self.deltaW>180, self.deltaW - 360, self.deltaW)
            self.sigma = np.where(self.sigma>180, self.sigma - 360, self.sigma)
            self.theta = np.where(self.theta>180, self.theta - 360, self.theta)

    def reOrderColumns(self):
        column_names = self.dataHsup.columns.tolist()  # Get the column names and move the last column to the second place
        last_column = column_names.pop(-1)             # Remove the last column from the list
        column_names.insert(1, last_column)            # Insert the last column in the second position
        self.dataHsup = self.dataHsup[column_names]    # Reorder the columns using the modified column_names list

    def print(self):
        print(self.dataHsup)

    def filterHdata(self, elimits = (0, 1), wlimits = (-180, 360), slimits = (-180, 360), surfaceNumber = 1):
        # Box limits to be preserved:
        emin = elimits[0]
        emax = elimits[1]
        wmin = wlimits[0]
        wmax = wlimits[1]
        smin = slimits[0]
        smax = slimits[1]
        if surfaceNumber == 1:
            self.dataHsup = self.dataHsup[ (((self.e1>=emin)&(self.e1<=emax))&((self.deltaW>=wmin)&(self.deltaW<=wmax))&((self.sigma>=smin)&(self.sigma<=smax))) ]
        if surfaceNumber == 2:
            self.dataHsup = self.dataHsup[ (((self.e2>=emin)&(self.e2<=emax))&((self.deltaW>=wmin)&(self.deltaW<=wmax))&((self.sigma>=smin)&(self.sigma<=smax))) ]

        #  Update variables from filtered dataframe:
        self.e1 = self.dataHsup['e1'].values
        self.e2 = self.dataHsup['e2'].values
        self.deltaW = self.dataHsup['deltavarpi'].values % 360
        self.sigma = self.dataHsup['sigma1'].values % 360
        self.H = -self.dataHsup['R'].values # Here it could be added the H_0 (Keplerian) term in the future
        self.theta = (self.sigma - (self.K2-self.K1)*self.deltaW) % 360
        # MAP MINDIS (OR DELTA) SURFACE?
        if self.mapMindis: self.mindis = self.dataHsup['mindis'].values
        # ANGLES DOMAIN ADEQUATION:
        self.anglesDomain(self.positiveAngles)
        self.reOrderColumns()

    def interpolateHdata(self):

        e_arr = np.unique(self.e1)
        #e_arr_aux = np.interp(np.arange(len(e_arr)), np.arange(len(np.unique(self.e2))), np.unique(self.e2))
        #e_arr_aux = interpolateArray(e_arr, np.unique(self.e2))
        e_arr_aux = np.linspace(np.min(np.unique(self.e2)), np.max(np.unique(self.e2)), len(e_arr))[::-1]
        w_arr = np.unique(self.deltaW)
       
       # PREPARO LOOP

        de = abs(e_arr[2] - e_arr[1])
        dw = abs(w_arr[2] - w_arr[1])

        Dnews = pd.DataFrame(columns = ['e1', 'e2', 'deltavarpi', 'sigma1', 'R', 'mindis'])

        # LOOP:
        for i_aux, e in enumerate(e_arr):
            e_aux = e_arr_aux[i_aux]

            for w in w_arr:

        # PUNTO ACTUAL
                Dsig = self.dataHsup[ (self.e1 == e) & (self.deltaW == w) ]
                sigmas = Dsig['sigma1'].values # % 360
                Hs = Dsig['R'].values
                minds = Dsig['mindis'].values
                # cantidad de sigmas:
                Nsig = len(sigmas)

        # HAY 3 INSTANCIAS DE INTERPOLACION POR PUNTO, HACIA X=e, HACIA Y=w Y HACIA AMBAS DIRECCIONES.

                DsigX = self.dataHsup[ (self.e1 == round(e + de, 3)) & (self.deltaW == round(w, 1)) ]
                if not DsigX.empty: 
                    sigmasX = DsigX['sigma1'].values % 360
                    HsX = DsigX['R'].values
                    mindsX = DsigX['mindis'].values
                    # cantidad de sigmas:
                    NsigX = len(sigmasX)

                    if (NsigX == Nsig):

                        for i in range(Nsig):

                            # punto actual
                            sigma = sigmas[i]
                            H     = Hs[i]
                            mind  = minds[i]

                            # punto siguiente
                            ind = find_nearest_ang(sigmasX, sigma) if (Nsig>1) else i 
                            sigmaX = sigmasX[ind]
                            HX = HsX[ind]
                            mindX = mindsX[ind]
                                
                            # punto nuevo
                            #sigma_new = ((sigma + sigmaX)/2) % 360
                            sigma_new = avg_ang(sigma,sigmaX)
                            if self.positiveAngles: sigma_new = np.where(sigma_new>180, sigma_new - 360, sigma_new)
                            H_new = (H+HX)/2
                            mind_new = (mind+mindX)/2

    #                        punto_new = [e + de/2, w, sigma_new, H_new, mind_new]
                            punto_new = [round(e + de/2, 3), round(e_aux + de/2, 3), round(w, 1), np.round(sigma_new, 1), H_new, mind_new]
                            Dnews.loc[len(Dnews)] = punto_new
            
                DsigY = self.dataHsup[ (self.e1 == round(e, 3)) & (self.deltaW == round(w + dw, 1)) ]   
                if not DsigY.empty: 
                    sigmasY = DsigY['sigma1'].values % 360
                    HsY = DsigY['R'].values
                    mindsY = DsigY['mindis'].values
                    # cantidad de sigmas:
                    NsigY = len(sigmasY)

                    if (NsigY == Nsig):

                        for i in range(Nsig):

                            # punto actual
                            sigma = sigmas[i]
                            H     = Hs[i]
                            mind  = minds[i]

                            # punto siguiente
                            #ind = find_nearest_ang(sigmasY, sigma) if (Nsig>1) else i 
                            ind = find_nearest_ang(sigmasY, sigma)
                            sigmaY = sigmasY[ind]
                            HY = HsY[ind]
                            mindY = mindsY[ind]

                            # punto nuevo
                            #sigma_new = ((sigma + sigmaY)/2) % 360
                            sigma_new = avg_ang(sigma,sigmaY)
                            if self.positiveAngles: sigma_new = np.where(sigma_new>180, sigma_new - 360, sigma_new)
                            H_new = (H+HY)/2
                            mind_new = (mind+mindY)/2

    #                         punto_new = [e, w + dw/2, sigma_new, H_new, mind_new]
                            punto_new = [round(e, 3), round(e_aux, 3), round(w + dw/2, 1), np.round(sigma_new, 1), H_new, mind_new]
                            Dnews.loc[len(Dnews)] = punto_new

                DsigXY = self.dataHsup[ (self.e1 == round(e + de, 3)) & (self.deltaW == round(w + dw, 1)) ]
                if not DsigXY.empty:     
                    sigmasXY = DsigXY['sigma1'].values % 360
                    HsXY = DsigXY['R'].values
                    mindsXY = DsigXY['mindis'].values
                    # cantidad de sigmas:
                    NsigXY = len(sigmasXY)

                    if (NsigXY == Nsig):

                        for i in range(Nsig):

                            # punto actual
                            sigma = sigmas[i]
                            H     = Hs[i]
                            mind  = minds[i]

                            # punto siguiente
                            #ind = find_nearest_ang(sigmasXY, sigma) if (Nsig>1) else i 
                            ind = find_nearest_ang(sigmasXY, sigma)
                            sigmaXY = sigmasXY[ind]
                            HXY = HsXY[ind]
                            mindXY = mindsXY[ind]

                            # punto nuevo
                            #sigma_new = ((sigma + sigmaXY)/2) % 360
                            sigma_new = avg_ang(sigma,sigmaXY)
                            if self.positiveAngles: sigma_new = np.where(sigma_new>180, sigma_new - 360, sigma_new)
                            H_new = (H+HXY)/2
                            mind_new = (mind+mindXY)/2

    #                         punto_new = [e + de/2, w + dw/2, sigma_new, H_new, mind_new]
                            punto_new = [round(e + de/2, 3), round(e_aux + de/2, 3), round(w + dw/2, 1), np.round(sigma_new, 1), H_new, mind_new]
                            Dnews.loc[len(Dnews)] = punto_new

        # Join both dataframes into one:
        self.dataHsup = pd.concat([self.dataHsup, Dnews], ignore_index=True)

        # Obtain variables from merged dataframe:
        self.e1 = self.dataHsup['e1'].values
        self.e2 = self.dataHsup['e2'].values
        self.deltaW = self.dataHsup['deltavarpi'].values % 360
        self.sigma = self.dataHsup['sigma1'].values % 360
        self.H = -self.dataHsup['R'].values # Here it could be added the H_0 (Keplerian) term in the future
        self.theta = (self.sigma - (self.K2-self.K1)*self.deltaW) % 360
        # MAP MINDIS (OR DELTA) SURFACE?
        if self.mapMindis: self.mindis = self.dataHsup['mindis'].values
        # ANGLES DOMAIN ADEQUATION:
        self.anglesDomain(self.positiveAngles)
        self.reOrderColumns()


