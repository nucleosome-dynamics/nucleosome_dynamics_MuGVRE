Nucleosome Dynamics Server  
------------  
  
This repository **wrapps Nucleosome Dynamics CLI** as a tool to make it accessible at the **MuGVRE platform** (http://multiscalegenomics.eu).  
  
#### Nucleosome Dynamics  
  
Nucleosome Dynamics is a suite of programs for defining nucleosome architecture and dynamics from MNase-seq and ATAC-seq experiments. They are implemented as a set of R libraries ([nucleR] (http://mmb.pcb.ub.es/gitlab/NuclDynamics/nucleR) and [nucDyn](http://mmb.pcb.ub.es/gitlab/NuclDynamics/NucleosomeDynamics_core)) and programs  ([Nucleosome Dynamics CLI](http://mmb.pcb.ub.es/gitlab/NuclDynamics/distpkg)) on top of which several implementations have been created in order to improve accessibility and usability. These include:  
  
1. MuGVRE platform: the domain-specific virtual research environment for 3D/4D genomics  - http://multiscalegenomics.eu  
2. Galaxy platform {[GIT](http://mmb.pcb.ub.es/gitlab/NuclDynamics/galaxy)}: the popular analysis platform - https://dev.usegalaxy.es  
3. Docker container {[GIT](http://mmb.pcb.ub.es/gitlab/NuclDynamics/docker)}: find the image [here](https://hub.docker.com/r/mmbirb/nucldyn)  
4. Singularity container {[GIT](https://github.com/nucleosome-dynamics/nucleosome_dynamics_singularity)}: find the image [here](https://www.singularity-hub.org/collections/2579)

#### MuGVRE platform

MuG Virtual Research Environment (MuGVRE) is the virtual research environment specifically developed for the [Multiscale complex Genomics (MuG)](http://multiscalegenomics.eu/MuG/) group, a research community focused on the 3D/4D genomics field. MuGVRE is a cloud-based e-infrastructure that integrates in a single framework rellevant data, tools and visualizers for the community in order to build a one-stop research platform.

> Learn more and access MuGVRE : http://multiscalegenomics.eu

'Nucleosome Dynamics' is one of the tools offered at the platform (check the [complete catalog](http://multiscalegenomics.eu/MuGVRE/tools-catalog/)).


## Requirements
* Nucleosome Dynamics:
    * Nucleosome Dynamics CLI
        -  NucDyn
        -  nucleR
* Python Modules:
    * mg-tool-api logger

#####  Nucleosome Dynamics CLI
Follow the instructions detailed at the [Nucleosome Dynamics CLI](http://mmb.pcb.ub.es/gitlab/NuclDynamics/distpkg) repository on how to install the set of R scripts and dependencies.

##### mg-tool-api
The present Nucleosome Dynamics wrapping makes use of some MuG-developed utilities like [mg-tool-api logger](https://github.com/Multiscale-Genomics/mg-tool-api/tree/master/utils), already included in this repository. It is a python logging falicity to report execution progress, exceptions, etc 


# Nucleosome Dynamics Server

This repository contains the python wrappers and the configuration files necessary to register and integrate 'Nucleosome Dynamics CLI' into  MuGVRE server:

- `nucleosome_dynamics.py` : MuGVRE tool. Interface for MuGVRE backend that receives web-user input files and calls 'Nucleosome Dynamics CLI' with them.
- `nucleosome_dynamics.json`: Tool definition. MuGVRE registry entry that defines 'Nucleosome Dynamics' tool with its input and output files, parameter definitions, resources, metadatam, etc.
- `test/`: Files for reproducing a 'Nucleosome Dynamics' run.

MuGVRE [documentation](http://multiscalegenomics.eu/MuGVRE/instructions/) includes a detailed explanation on the number of steps necessary to integrate a new tool into the platform.

##### Running 'Nucleosome Dynamics' tool

'Nucleosome Dynamics' tool is mean to be executed though the MuGVRE web interface. Yet, for checking this code integrity, the 'test/' folder contains a simple bash script and the required sample data that emulates how  MuGVRE server deals with a tool execution. 

MUGVRE platform interacts with a tool through a main wrapper, here `nucleosome_dynamics.py`. It always take as inputs a set of configuration JSON files that contains the rellevant information gathered from the web-user (input files and arguments) for a particular run. Also includes the data for the expected output files. For nucleosome dynamics, the script sequencially runs 'Nucleosome Dynamics CLI' according the given configuration files to end up executing a tailored Nucleosome Dynamics workflow based on the web-user indications.

Arguments
* config: JSON file containing the run parameters, working directory, analyses to be performed, and their specific arguments 
* in_metadata: JSON file containing MuG metadata attached to each of the web-user input files
* out_metadata: Filename for the JSON that will contain output file metadata
* log_file: Filename for the log file

The following bash script runs the MuGVRE tool (nucleosome_dynamics.py) with the the configuration files set to read the MNase-seq sample data found at `test/data`.

```sh
cd test/test
bash test_0_AnalyseMNaseseqdata.sh
```

The `test/test` content has been automatically created by MuGVRE during the tool integration procedure (http://multiscalegenomics.eu/MuGVRE/instructions/) in order to check that VRE prepares the right configuration files based on the given tool definition. In this way, tool developers can easily check whether the MUGVRE tool (here nucleosome_dynamics.py) correctly fits the MuGVRE job processing.

