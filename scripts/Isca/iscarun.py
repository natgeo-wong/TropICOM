import os
import numpy as np
import f90nml
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE, GFDL_DATA

NCORES = 64
RESOLUTION = 'T85', 96

base_dir = os.path.dirname(os.path.realpath(__file__))
cb = IscaCodeBase.from_directory(GFDL_BASE)

cb.compile()

exp = Experiment("IscaZonal", codebase=cb)

diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'zsurf')
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'height_half', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('atmosphere', 'precipitation', time_avg=True)

exp.clear_rundir()
exp.diag_table = diag
exp.namelist = f90nml.read('namelist.nml')
exp.set_resolution(*RESOLUTION)

if __name__=="__main__":
    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2,21):
        exp.run(i, num_cores=NCORES)