# nucleServ

Scripts and wrappers to run nucleR and NucleosomeDynamics on the backend of a server

IMPORTANT: note that the paths to the sourced R functions are hardcoded


usage: pipeline.py [-h] -c config_file [-a] [--preproc1] [--preproc2]
                   [--nucleR1] [--nucleR2] [--nucdyn]

Wrapper to run nucleR and NucleosomeDynamics

optional arguments:
  -h, --help            show this help message and exit
  -c config_file, --config config_file
                        configuration file
  -a, --all             run the whole pipeline

single calculations:
  --preproc1            preprocessing of the first experiment
  --preproc2            preprocessing of the second experiment
  --nucleR1             nucleR on the first experiment
  --nucleR2             nucleR on the second experiment
  --nucdyn              NucleosomeDynamics
