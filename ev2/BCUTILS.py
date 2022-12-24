import matplotlib.pyplot as plt; import numpy as np ;import pandas as pd
import pwlf ;import h5py
#Read data from web API and monitor HTTP traffic 
import requests  
#specialised functions to query Illustris TNG data 
import illustris_python as il
#specific functions for fitting utilities
from scipy.optimize import curve_fit ; from scipy.signal import medfilt, savgol_filter
#
#Important constants
#

headers = {"api-key":"849c96a5d296f005653a9ff80f8e259e"}

#
#Main Code
#

#----------------------------------------------------------------------------------------------------------------------------------------
#CLASS UTILITY
#----------------------------------------------------------------------------------------------------------------------------------------
class UTILITY:
    def __init__(self) -> None:
        pass
    
    def get(path, params = None):
        r'''
        func get(path, params = None)
        
        Utility function to read data from API using http requests: custom modification expanding upon requests.get():
        
        
        Expansions:
        
        Content Filtering: built to read information from illustris TNG API:
        Error Raising: include raise_for_status function to display all error codes (expect HTTP return 200 to indicate all OK)
        
        Valid data types =
        
        ->application/json
        
        -> .hdf5 
        
        
        '''
        #utility function for API reading 

        #Make API request - 
        # Path: url to api page 
        #Params - misc ; Headers = api key 
        r = requests.get(path, params=params, headers=headers)

        #HTTP code - raise error if code return is not 200 (success)
        r.raise_for_status()
        
        #detect content type (json or hdf5) - run appropriate download programme
        
        if r.headers['content-type'] == 'application/json':
            return r.json() # parse json responses automatically

        if 'content-disposition' in r.headers:
            filename = r.headers['content-disposition'].split("filename=")[1]
            with open(filename, 'wb') as f:
                f.write(r.content)
            return filename # return the filename string

        return r

    def line(m,x,b):
        '''
        straight line function y=m*x+b
        '''
        y = 10**((m*x)+b)
        return y 
    
    def linear_fit(a,x,b):
        f = (a*x)+b
        return f

    def sq_fit(x,a,b,c):
        '''
        quadratic function (y=ax**2+b*x+c)
        '''
        f = (a*(x**2))+(b*x)+c
        return f

#----------------------------------------------------------------------------------------------------------------------------------------
#CLASS SNAPTOOLS
#----------------------------------------------------------------------------------------------------------------------------------------
class snaptools:
    def __init__(self) -> None:
        pass
    
    def classification(snapID,simID):
        baseurl = "https://www.tng-project.org/api/{}/snapshots/{}/subhalos/".format(simID,snapID)
        inq = "?sfr__gt=0.0"
        
        counturl = baseurl + inq
        getcount = UTILITY.get(counturl)
        count = getcount['count']
        
        nextq = "?limit={}&sfr__gt=0.0".format(count)
        suburl = baseurl+nextq
        snapshot = UTILITY.get(suburl)
        subIDS = [snapshot['results'][i]['id'] for i in range(snapshot['count'])]
        mass = []; sfr=[]
        for i,id in enumerate(subIDS):
            mass.append(snapshot['results'][i]['mass_log_msun'])
            sfr.append(snapshot['results'][i]['sfr'])

        df = pd.DataFrame({
            "id": subIDS,
            "mass": mass,
            "sfr": sfr
        })
        fpath = "starforming_subhalos_snapshot{}.csv".format(snapID)
        df.to_csv(fpath)
        return print("Subhalos classified for snapshot {}".format(snapID))
    
    def datafilters(df,quant,upper,lower):
        df = df.copy()
        df = df[df[quant]<upper]
        df = df[df[quant]>lower]
        return df

#----------------------------------------------------------------------------------------------------------------------------------------
#CLASS PLOTS
#----------------------------------------------------------------------------------------------------------------------------------------
class plots:
    def __init__(self,subID):
        self.subID = subID
    
    def metgrad(self,dfin):
        df = dfin.sort_values(by='rad').copy()
        plt.figure(figsize=(20,12))
        plt.plot(df['rad'],12+np.log10(df['met']),'g+')
        plt.xlabel("Radius (Kpc)") ; plt.ylabel("12+$log_{10} ({O}/{H}$")
        fpath = "metparticles_{}.png".format(self.subID)
        plt.savefig(fpath)
        plt.close()
        
    def metgradsfr(self,dfin):
        df = dfin.sort_values(by='rad').copy()
        plt.figure(figsize=(20,12))
        plt.scatter(df['rad'],12+np.log10(df['met']),c=df['sfr'],cmap = 'viridis')
        plt.xlabel("Radius (Kpc)") ; plt.ylabel("12+$log_{10} ({O}/{H}$")
        fpath = "metparticlessfr_{}.png".format(self.subID)
        plt.savefig(fpath)
        plt.close()
    
    def linearmetgrad(self,dfin):
        df = dfin.sort_values(by='rad').copy()
        x0 = list(df['rad'])
        interp = savgol_filter(12+np.log10(df['met']),window_length=21,polyorder=5)
        popt,pcov = curve_fit(UTILITY.linear_fit, df['rad'],(12+np.log10(df['met'])),sigma = 1/df['sfr'],absolute_sigma=True)
        
        plt.figure(figsize=(20,12))
        plt.plot(df['rad'], interp, 'r-')
        plt.plot(df['rad'], UTILITY.linear_fit(x0,*popt),'g--')
        plt.xlabel("Radial position (Kpc)")
        plt.ylabel("12+$log_{10} ({O}/{H}$")
        fpath = "linearfit_{}.png".format(self.subID)
        plt.savefig(fpath)
        plt.close()
        
    def singlebreak(self,dfin,breakpoint):
        df = dfin.copy()
        df.sort_values(by="rad",inplace = True)
        interp = savgol_filter(12+np.log10(df['met']),window_length=21,polyorder=5)
        x0 = np.array([min(df['rad']), breakpoint, max(df['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met2']),weights=1/df['sfr'])
        my_pwlf.fit_with_breaks(x0)
        xHat = np.linspace(min(df['rad']), max(df['rad']), num=10000); yHat = my_pwlf.predict(xHat)
        
        plt.figure(figsize=(20,12))
        plt.plot(df['rad'], interp, 'b--')
        plt.plot(xHat,yHat, 'g-')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        fpath = 'singlebreak_{}'.format(self.subID)
        plt.savefig(fpath)
        plt.close()

    def doublebreak(self,dfin,break1,break2):
        df = dfin.copy()
        df.sort_values(by="rad",inplace = True)
        interp = savgol_filter(12+np.log10(df['met']),window_length=21,polyorder=5)
        x0 = np.array([min(df['rad']), break1,break2, max(df['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met']),weights=1/df['sfr'])
        my_pwlf.fit_with_breaks(x0)
        xHat = np.linspace(min(df['rad']), max(df['rad']), num=10000); yHat = my_pwlf.predict(xHat)
        
        plt.figure(figsize=(20,12))
        plt.plot(df['rad'], interp, 'b--')
        plt.plot(xHat,yHat, 'g-')
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        fpath = 'doublebreak_{}.png'.format(self.subID)
        plt.savefig(fpath)
        plt.close()
        
    def gas_visual(self,dfin,decp):
                df_valid = dfin.round(decp)
                df_valid = df_valid.groupby(['x','y'])['dens'].sum().reset_index()
                plt.figure(figsize=(20,12), dpi=500)
                plt.style.use('dark_background')
                plt.scatter(df_valid['x'],df_valid['y'],c=(np.log10(df_valid['dens'])),cmap='inferno',
                            vmin=(min(np.log10(df_valid['dens']))),vmax =(0.7*max(np.log10(df_valid['dens']))))
                plt.xlabel('$\Delta x$ [kpc/h]')
                plt.ylabel('$\Delta y$ [kpc/h]')
                plt.colorbar(label='log10(Gas Mass)')
                filename = 'Mgass_sub_{}.png'.format(self.subID)
                plt.savefig(filename)
                plt.close()
        
    def met_histogram(self,dfin,extra):
        df = dfin.copy()
        df.sort_values(by='rad',inplace=True)
        df.to_csv("historgram.csv")
        plt.figure(figsize=(20,12))
        plt.hist2d(df['rad'],12+np.log10(df['met']),bins=[200,200], weights=1/df['sfr'],cmap='PuOr')
        if (extra=='Y'):
            interp = savgol_filter(12+np.log10(df['met']),window_length=21,polyorder=5)
            plt.plot(df['rad'],interp,'b-')
        else:
            print("No Sav Filter Apploied")
        plt.xlabel("Radius (Normalised Code Units)")
        plt.ylabel("12+$log_{10}$ $(O/H)$")
        filename = 'historgram_met_{}.png'.format(self.subID)
        plt.savefig(filename)
        plt.close()

#----------------------------------------------------------------------------------------------------------------------------------------
# CLASS BOOTSTRAPPING
#----------------------------------------------------------------------------------------------------------------------------------------     
class bootstrapping:
    def __init__(self,subID):
        self.subID = subID
        
    def linear(self,dfin,runs,frac):
        avals = []; bvals = [] ; pcov0 = []; pcov1=[]
        df = dfin.copy()
        for i in range(runs):
            sample = df.sample(frac=frac, replace=False)
            popt,pcov = curve_fit(UTILITY.linear_fit,sample['rad'],sample['met'])
            avals.append(popt[0])
            bvals.append(popt[1])
            pcov0.append(pcov[0])
            pcov1.append(pcov[1])
        return avals,bvals,pcov0, pcov1
    
    def singlebreak(self,dfin,runs,frac,break1):
        df = dfin.copy()
        inner = []; outer= []
        for i in range(runs):
            df = df.sample(frac=frac, replace=False)
            x0 = np.array([min(df['rad']), break1, max(df['rad'])])
            my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met2']),weights=1/df['sfr'])
            my_pwlf.fit_with_breaks(x0)
            inner.append(my_pwlf.slopes[0])
            outer.append(my_pwlf.slopes[1])
        
        return inner,outer
    
    def doublebreak(self,dfin,runs,frac,break1,break2):
            df = dfin.copy()
            l1=[];l2=[];l3=[]
            for i in range(runs):
                df = df.sample(frac=frac, replace=False)
                x0 = np.array([min(df['rad']), break1,break2, max(df['rad'])])
                my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met']),weights=1/df['sfr'])
                my_pwlf.fit_with_breaks(x0)
                slope1 = my_pwlf.slopes[0]
                slope2 = my_pwlf.slopes[1]
                slope3 = my_pwlf.slopes[2]
                l1.append(slope1);l2.append(slope2);l3.append(slope3)
            
            return(slope1,slope2,slope3)


#----------------------------------------------------------------------------------------------------------------------------------------
# CLASS CSVTOOLS
#----------------------------------------------------------------------------------------------------------------------------------------
class csvtools:
    def __init__(self) -> None:
        pass
    
    def docombine(i):
        try:
            path1 = "files/historycutouts/evdir_{}/slope{}.csv".format(i,i)
            path2 = "files/historycutouts/evdir_{}/historydata_{}.csv".format(i,i)
            df1 = pd.read_csv(path1)
            df2 = pd.read_csv(path2)
            df2.insert(6,"slope1",df1['slope1'])
            df2.insert(7,"slope2",df1['slope2'])
            df2.insert(8,"slope3",df1['slope3'])
            path = "files/historycutouts/evdir_{}/alldata{}.csv".format(i,i)
            df2.to_csv(path)
            print("done for {}".format(i))
        except FileNotFoundError as e:
            print("no file for {}".format(i))
            return print(e)
#----------------------------------------------------------------------------------------------------------------------------------------
#CLASS TREE    
#----------------------------------------------------------------------------------------------------------------------------------------
class tree:
    def __init__(self,subhalo):
        self.ID = subhalo
        baseurl = "https://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/"
        self.suburl = baseurl + str(self.ID)
        print("done for subhalo {}".format(subhalo))
        
    def treeget(path,dir, params = None):
        #utility function for API reading 

        #Make API request - 
        # Path: url to api page 
        #Params - misc ; Headers = api key 
        r = requests.get(path, params=params, headers=headers)

        #HTTP code - raise error if code return is not 200 (success)
        r.raise_for_status()
        
        #detect content type (json or hdf5) - run appropriate download programme
        
        if r.headers['content-type'] == 'application/json':
            return r.json() # parse json responses automatically

        if 'content-disposition' in r.headers:
            filename = str(dir)+r.headers['content-disposition'].split("filename=")[1]
            with open(filename, 'wb') as f:
                f.write(r.content)
            return filename # return the filename string

        return r
    
    def get_tree(self,dir):
        sub = UTILITY.get(self.suburl)
        self.mpb = self.treeget(sub['trees']['sublink_mpb'],dir)
        return print("Tree got and saved to {}".format(str(self.dir)))

    def history_trace(self):
        subID = self.ID
        mpb = self.mpb
        with h5py.File (mpb,'r') as f:
            self.snapshots = list(f['SnapNum'][:])
            self.subhalos= list (f['SubfindID'][:])
        
        self.snapshots.reverse(); self.subhalos.reverse()
        df = pd.DataFrame({"snapshot":self.snapshots,"subhalo":self.subhalos})
        
    def snapshot_filter(self,dfin,filterlist):
        df=df.copy()
        df=df[df['snapshot'].isin(filterlist)]
        return df
    
#----------------------------------------------------------------------------------------------------------------------------------------
# CLASS SHAPEANALYSIS
#----------------------------------------------------------------------------------------------------------------------------------------
class shapeanalysis:
    def __init__(self) -> None:
        pass
    
    def get_char_3D(x,y,z):
        '''
        get_char_3d(x,y,z)
        
        INPUTS:
        
        x: float | Inner gradient of subhalo metallicity gradient
        
        y: float | Central gradient of subhalo metallicity gradient
        
        z: float | Outer gradient of subhalo metallicity gradient 
        
        
        Outputs: int range(1-7) | identifying number describing gradient broken fit shape according to 
                ->shape 1:  steep negative inner followed by shallow negative mid/outer
                ->shape 2:  shallow negative inner followed by steeper negative mid/outer
                ->shape 3:  positive inner followed by negative outers
                ->shape 4:  positive inner followed by positive outer( Something has gone wrong)
                ->shape 5:  steep inner down followed by shallow negative outer
                ->shape 6:  steep inner up followed by shallow negative outer
                ->shape 7:  No shape criteria met

        '''
        
        
        #shape 1: steep negative inner followed by shallow negative mid/outer
        if x<0 and y<0 and abs(x)<abs(y):
            return 1
        #shape 2: shallow negative inner followed by steeper negative mid/outer
        elif x<0 and y<0 and abs(y)<abs(x):
            return 2
        #shape 4:positive inner followed by positive outer( Something has gone wrong)
        elif x>0 and y>0:
            return 4
        #shape 5: steep inner down followed by shallow negative outer:
        elif x<-0.5 and y<0 and y>-0.1:
            return 5
        #shape6: steep inner up followed by shallow negative outer:
        elif x>0.5 and y<0 and y>-0.1:
            return 6
        #shape 3: positive inner followed by negative outers:
        elif x>0 and y<0:
            return 3
        else:
            return 7
        
    def get_char_2D(x,y):
        #shape 1: steep negative inner:
        if x<-0.1 and y<0:
            return 1 
        #shape 2: steep positive inner:
        elif x>0.1 and y<0:
            return 2
        #shape 3: anything else 
        else:
            return 3 
    
    def get_inner_direction(x):
        '''
        Inputs: 
        
        X: float (num)
        
        Gradient of inner slope of metallicity gradient fit 
        
        outputs: Shape characteristic identifying number 
        
        '''
        if x>0:
            return 1
        elif x<0: 
            return 2
        else:
            return 3    
        
#----------------------------------------------------------------------------------------------------------------------------------------
# CLASS CUTOUT_SUBHALO
#----------------------------------------------------------------------------------------------------------------------------------------
class cutout_subhalo:
    def __init__(self,subID,snapID,simID,primeID):
        
        self.subID = subID
        self.snapID = snapID
        self.simID = simID
        self.primeID = primeID
        redurl = 'http://www.tng-project.org/api/'+str(simID)+'/snapshots/'+str(snapID)
        utility = UTILITY.get(redurl)
        redshift = utility['redshift']
        self.redshift =redshift
        hubble = 0.7
        scalefac = 1./(1.+redshift) #calculate scale factor

        #
        #investigate particle level subhalo data 
        #
        # begin by downloading subhalo cutout       
        self.subURL = 'http://www.tng-project.org/api/TNG50-1/snapshots/{}/subhalos/{}/'.format(str(snapID),str(subID))
        subhalo = UTILITY.get(self.subURL)
        #
        #obtain global properties
        #
        self.mass = subhalo['mass_log_msun']
        self.tot_sfr = subhalo['sfr']
        self.tot_met = subhalo['gasmetallicity']
        self.Rhalf = subhalo['halfmassrad']
        self.crit_dist = 5* self.Rhalf
        self.stellarphotometrics = subhalo['stellarphotometricsrad']
        
        
        self.centre = np.array([subhalo['pos_x'],subhalo['pos_y'],subhalo['pos_z']])
        
        cutout ="files/historycutouts/evdir_{}/cutout_{}.hdf5".format(self.primeID,self.subID)
        with h5py.File(cutout,'r') as f:
            sfr = f['PartType0']['StarFormationRate'][:]
            co_ords = f['PartType0']['Coordinates'][:]
            hcoldgas  = np.where( (sfr > 0.0))[0]
            self.pgas_coo = f['PartType0']['Coordinates'][hcoldgas]
            self.pgas_m = f['PartType0']['Masses'][hcoldgas]
            self.pgas_vel = f['PartType0']['Velocities'][hcoldgas]
            self.pgas_met = f['PartType0']['GFM_Metallicity'][hcoldgas]
            self.pgas_sfr = f['PartType0']['StarFormationRate'][hcoldgas]
        self.test = len(hcoldgas)

        self.pgas_coo -= self.centre[None,:]

    def align_dfgen(self):
        '''
        Inputs: Self | Particle level data input from object creation
        
        Computes angular momentum of subhalo gas, uses information to create a transformation matrix to align the subhalo such that x,y plane is
        perpendicular to z 
        
        '''
        _coo = np.copy(self.pgas_coo)
        _vel = np.copy(self.pgas_vel)
        _m = np.copy(self.pgas_m)
        
        self.ang_mom_3D = np.sum(_m[:,None,]*np.cross(_coo,_vel),axis=0)
        self.ang_mom = self.ang_mom_3D/ np.sum(_m)

        j=self.ang_mom/np.linalg.norm(self.ang_mom)
        #normalised specific angular momentum 
        
        x = np.array([1,2,3])
        x = x-(x.dot(j)*j) #make x orthogonal to j
        
        x/= np.linalg.norm(x) # normalise
        
        y = np.cross(j,x)#create 3rd vector - orth to x,j
        
        
        A = (x,y,j) # transformation matrix
        
        self.pgas_coo=np.dot(A,self.pgas_coo.T).T # change co-ordinates
        self.pgas_vel = np.dot(A,self.pgas_vel.T).T
        self.radial = np.sqrt((self.pgas_coo[:,0]**2)+(self.pgas_coo[:,1]**2))
        
        df = pd.DataFrame({
            "x":self.pgas_coo[:,0],
            "y":self.pgas_coo[:,1],
            "z":self.pgas_coo[:,2],
            "rad": self.radial,
            "mass":self.pgas_m,
            "met":(self.pgas_met),
            "sfr":self.pgas_sfr
            })
        df=df[df['rad']<5*self.Rhalf]
        f=df.sample(frac=0.1,replace=False)

        #print(df)
        self.df = df
        
    def filter(self):
        df = self.df.copy()
        spr = self.stellarphotometrics
        z_max = 0.1*spr
        df = df[df['z']<z_max]
        df.rad =10*((df.rad-df.rad.min())/(df.rad.max()-df.rad.min()))
        self.df = df
        
        return df
    
    def linearfit(self, dfin):
        dfin.sort_values(by='rad',inplace=True)
        popt,pcov = curve_fit(UTILITY.linear_fit, dfin['rad'],np.log10(dfin['met'])+12,sigma=1/dfin['sfr'],absolute_sigma=True)
        
        return popt[0]
    
    def piecewise(self,dfin,breakpoint):
        df = dfin.copy()
        df.sort_values(by="rad",inplace = True)
        x0 = np.array([min(df['rad']), breakpoint, max(df['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met']),weights=1/df['sfr'])
        my_pwlf.fit_with_breaks(x0)
        slope1 = my_pwlf.slopes[0]
        slope2 = my_pwlf.slopes[1]
        
        return (slope1,slope2)
    
    def doublepiecewise(self,dfin,breakpoint1,breakpoint2):
        df = dfin.copy()
        df.sort_values(by="rad",inplace = True)
        x0 = np.array([min(df['rad']), breakpoint1,breakpoint2, max(df['rad'])])
        my_pwlf = pwlf.PiecewiseLinFit(df['rad'], 12+np.log10(df['met']),weights=1/df['sfr'])
        my_pwlf.fit_with_breaks(x0)
        slope1 = my_pwlf.slopes[0]
        slope2 = my_pwlf.slopes[1]
        slope3 = my_pwlf.slopes[2]
        
        return (slope1,slope2,slope3)
    
#----------------------------------------------------------------------------------------------------------------------------------------
# CLASS DODIRECTORY
#----------------------------------------------------------------------------------------------------------------------------------------
class dodirectory:
    def __init__(self,primeID):
        self.primeID = primeID
    
    def getlist(self,fpath):
        df = pd.read_csv(fpath)
        snapshots = list(df['snapshots'])
        subhalos = list(df['subhalos']) 
        return snapshots,subhalos
    
#----------------------------------------------------------------------------------------------------------------------------------------
# CLASS CUTOUT_PROCESSING
#----------------------------------------------------------------------------------------------------------------------------------------
class cutout_processing:
    def dosingle(sub,snap,prime):
        subhalo = cutout_subhalo(sub,snap,'TNG50-1', prime)
        subhalo.align_dfgen()
        df= subhalo.filter()
        subID = sub
        snapID = snap
        slope1,slope2,slope3= subhalo.doublepiecewise(df,3,8)
        print("done for subhalo {} snapshot {}".format(sub,snap))
        return (subID, snapID,slope1,slope2,slope3)
    
    def dodir(i):
        try:
            data = dodirectory(i)
            snapshots,subhalos = data.getlist()
            subs = [];snaps=[];s1=[];s2=[];s3=[]
            for j in range(len(snapshots)):
                subID,snapID,slope1,slope2,slope3 = dosingle(subhalos[j],snapshots[j],i)
                subs.append(subID);snaps.append(snapID)
                s1.append(slope1);s2.append(slope2);s3.append(slope3)
            df = pd.DataFrame({
                'subhalo':subs,
                'snapshot':snaps,
                'slope1':s1,
                'slope2':s2,
                'slope3':s3  
            })
            fpath = "files/historycutouts/evdir_{}/slope{}.csv".format(i,i)
            df.to_csv(fpath)
            return print("done for descendant {}".format(i))
        except OSError as e:
            return print(e)
        except TypeError as e:
            return print(e)
        except KeyError as e:
            return print(e)
        except IndexError as e:
            return print(e)
        except ValueError as e:
            return print(e)
        
