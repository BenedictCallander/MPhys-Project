import h5py 
import glob
basePath = '/home/AstroPhysics-Shared/DATA/IllustrisTNG/TNG50-1/output/snapdir_099/*.hdf5'
fail=[]
for file in glob.glob(basePath):
    try : 
        test = [el for el in h5py.File(file,"r")]
    except OSError:
        fail.append(file)
print(fail)