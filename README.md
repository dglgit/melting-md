# melting

This is an attempt at simulating Differential Screening Fluorimetry using molecular dynamics. I used OpenMM to simulate and Ambertools to set up the simulation. Both are free and open source tools. 

`gen_cplx.tleap` is the config file for `tleap` to generation prmtop and inpcrd files, which essentially describe the system being simulated. Run it with `tleap -f gen_cplx.tleap`
`antechamber` is needed to generate the `.frcmod` file from the `.mol` file describing the ligand

Input simulation parameters in a json format to `meltv2.py`. See `sample.json` to see what parameters are needed.
Use the same json file as input to `trajAnal.py` to plot and analyze your data. 
