# scRNA-seq-pipeline
 scRNA-seq snakemake repository

# Dependencies
Dependencies (w/ versions) are listed in the environment.yml file. Install these dependencies using a conda environment. 

'''
conda env create -f /path/to/gitclone/environment.yml
'''

You will then want to activate the environment. 

'''
conda activate Pipeline_conda
'''

Once active you will also need to install the following R packages via an attached R script. 

'''
R R-packages-to-download.R
'''

Now you are almost ready to run our snakemake pipeline. Please download our data from here [_____] (google.com).

After downloading the data you will need to change the configuration file. 

'''
OUTDIR: /path/to/where/outputs/should/go #where output data from pipeline run should be stored
SNAKEDIR: /path/to/gitclone/scRNA-seq-pipeline #path to the github repository clone 
DATADIR: /path/to/where/downloaded/data/is #where the data you just downloaded is stored
'''






