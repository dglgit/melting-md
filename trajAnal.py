import mdtraj 
import numpy as np
import matplotlib.pyplot as plt
trajectory=mdtraj.load_dcd('melt_prot_305-325-1_traj.dcd',top='melt_prot_305-325-1_eqed.pdb')
#trajectory=mdtraj.load_dcd('melt_prot_310-360_traj.dcd',top='melt_prot_310-360_eqed.pdb')
#trajectory=mdtraj.load_dcd('melt_complex2_330-350_traj.dcd',top='melt_complex2_330-350_eqed.pdb')
#trajectory=mdtraj.load_dcd('melt_complex2_traj.dcd',top='melt_complex2_eqed.pdb')
rmsds=mdtraj.rmsd(trajectory,trajectory,0)
def movingAverage(x,window=10):
    res=np.empty_like(x)
    for i in range(len(x)):
        res[i]=x[max(0,i-window//2):min(len(x)-1,i+window//2)].mean()
    return res
#rmsds=movingAverage(rmsds)
temp=310
means=[]
temps=[]
stds=[]
for i in range(0,len(rmsds),50):
    roi=rmsds[i:i+50]
    mean=roi.mean()
    std=np.std(roi)
    print(f'temp: {temp}, mean: {mean}, std:{std}')
    
    means.append(mean)
    temps.append(temp)
    stds.append(std)
    temp+=1
#print(rmsds)
plt.plot(rmsds)
plt.show()
plt.scatter(temps,means)
plt.errorbar(temps,means,yerr=np.array(stds))
plt.show()
