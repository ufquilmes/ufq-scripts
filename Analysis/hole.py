import numpy as np
import MDAnalysis as mda
from mdahole2.analysis.hole import HoleAnalysis

pdb_file  = "../../GBSA/apt_1/solucion.prmtop"
traj_file = "md_complejo_1.nc"
hole_exec = "/home/aormazabal/miniconda3/envs/mdanalysis/bin/hole"
out_vmd   = "hole_1.vmd"

u = mda.Universe(pdb_file, traj_file)

ha = HoleAnalysis(u, executable=hole_exec)
ha.run()

min_r = ha.min_radius()
np.savetxt("min_radius.txt", min_r)

ha.create_vmd_surface(filename=out_vmd)

try:
    ha.plot()
except Exception:
    pass

