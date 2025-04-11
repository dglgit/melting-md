import mdtraj 
import numpy as np
import matplotlib.pyplot as plt
import json
import sys

if len(sys.argv)<2:
    print('provide config file')
    exit()
fname=sys.argv[1]
f=open(fname,'r')
d=json.load(f)
f.close()
base=d['output_root']
dcd=base+'_traj.dcd'
ref=base+'_eqed.pdb'
dcd_interval=d['dcd_interval']
eq_frames=d['eq_steps']//dcd_interval
eq_temp=d['eq_temperature']
temp_frames=d['steps_per_temp']//dcd_interval

trajectory=mdtraj.load_dcd(dcd,top=ref)
rmsds=mdtraj.rmsd(trajectory,trajectory,0)
rogs=mdtraj.compute_rg(trajectory)
def movingAverage(x,window=10):
    res=np.empty_like(x)
    for i in range(len(x)):
        res[i]=x[max(0,i-window//2):min(len(x)-1,i+window//2)].mean()
    return res
#rmsds=movingAverage(rmsds)

if d['start'] is not None:
    temps=list(np.arange(d['start'],d['end'],d['step']))
else:
    temps=d['temperatures']
rmsd_means=[rmsds[:eq_frames].mean()]
rog_means=[rogs[:eq_frames].mean()]
tx=[eq_temp]
rmsd_stds=[np.std(rmsds[:eq_frames])]
rog_stds=[np.std(rogs[:eq_frames])]
temps=[eq_temp]+temps
temp_idx=1
for i in range(eq_frames,len(rmsds),temp_frames):
    rmsd_roi=rmsds[i:i+temp_frames]
    rog_roi=rogs[i:i+temp_frames]
    rmsd_mean=rmsd_roi.mean()
    rog_mean=rog_roi.mean()
    rmsd_std=np.std(rmsd_roi)
    rog_std=np.std(rog_roi)
    #print(f'temp: {temp}, mean: {mean}, std:{std}')
    
    rmsd_means.append(rmsd_mean)
    rog_means.append(rog_mean)
    tx.append(temps[temp_idx])
    rmsd_stds.append(rmsd_std)
    rog_stds.append(rog_std)
    temp_idx+=1
#print(rmsds)
plt.plot(rmsds)
plt.title(f'rmsd over all frames. equilibration accounts for the first {eq_frames} frames')
plt.show()
plt.scatter(temps,rmsd_means)
plt.title("average RMSD over a temperature with 2 SD error bars")
plt.errorbar(temps,rmsd_means,yerr=np.array(rmsd_stds))
plt.show()
plt.scatter(temps,rog_means)
plt.errorbar(temps,rog_means,yerr=np.array(rog_stds))
plt.title('Radius of Gyration')
plt.show()
