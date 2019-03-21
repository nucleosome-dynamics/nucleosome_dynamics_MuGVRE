# nucleServ - Nucleosome Dynamics Server

Scripts and wrappers to run nucleR and NucDyn on the backend of a server

## Analyses

#### NucDyn

Primary data
* Position: region where a nucleosome movement is detected
* Type: change in the nucleosome map
* Score: magnitude of the change

Attributes
* class: type of hotspot (see help for all possible types)
* nreads: number of reads involved in this movement


#### nucleR

Primary data
* Score: Positionning score. It is calculated as the weighted sum of width and height scores.

Attributes
* score_width: Witdth score. It is a measure of how sharp a peak is. A value of 0 would be an extremely wide peak and a value of 1 a very sharp one.
* score_height: Height score. Tells how large a peak of a nucleosome is. The bigger this number, the higher the peak.
* class: Whether the nucleosome is well-positioned (W) or fuzzy (F) or undetermined. The taken value depends on score_h and score_w. Undetermined means the exact position of the nucleosome cannot be determined due to strong fuzziness.


#### Stiffness

Primary data
* Score: Stiffness estimation. It represents the energy required to move the nucleosome (expressed in kcal/mol/bp), it's derived from the Gaussian standard deviation.

Attributes
* nucleR_score: the nucleR score given to that nucleosome
* nucleR.class: the nucleR class given to that nucleosome
* gauss_k: the height of the peak of the gaussian curve
* gauss_m: the position of the peak of the gaussian curve
* gauss_sd: the standard deviation of the gaussian curve


#### Periodicity

Primary data
* Position: gene position (from TSS to TTS) 

Attributes
* nucleosome_first: First nucleosome of the gene.
* nucleosme_last: Last nucleosome of the gene.
* score_phase: Is a measure of the phase between the first and the last nucleosome. A score of 0 means the nucleosome are completely phased and a score of 82 corresponds to totally antiphased nucleosomes.
* score_autocorrelation: It is directly computed from the experimental coverage and is quantitative measure of the periodicity of nucleosomes inside the gene body.


#### TSS classes

Primary data
* Position: Region between the dyads of two nucleosomes surrounding the TSS.

Attributes
* classification: Descriptor of the Transcription Start Site. See the help for possible options.
* distance: Distance in base pairs between the nucleosome +1 and the nucleosome -1.
* nucleosome minus1: Position of the nucleosome -1.
* nucleosome plus1: Position of the nucleosome +1
* TTS_position: Position of the Transcription Start Site.

# Complete workflow

#### nucleosome_dynamics.py

Sequencially runs the analyses above listed in a complete Nucleosome Dynamics workflow.

Arguments
* config: JSON file containing workflow parameters
* in_metadata: JSON file containing MuG metadata files
* out_metadata: Filename for the JSON that will contain output file metadata
* log_file: Filename for the log file

The following bash script runs `nucleosome_dynamics.py` with the sample MNase-seq data found at `test/data`. 
```sh
cd test/
bash test_0_AnalyseMNaseseqdata.sh 
```

