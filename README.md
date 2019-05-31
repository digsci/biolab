# diatom-pipeline
Pipeline for diatom analysis

## Running Pipeline

Ensure you have Docker installed.  Instructions on how to do that are found at https://www.docker.com/

If you have built this container before please run `docker-compose build` before running `docker-compose up`


### If image has not been built before

1. In the directory containing these scripts, create a folder called **"sequences"** and then save all the fastq files in this folder.
2. Open a terminal in the directory containing all the scripts and run `docker-compose up`

**Make sure that the docker-compose.yml file is in the same directory as the other scripts and save all the fastq files you wish to analyse to the 'sequences' folder on your local machine / VM.**

### Once the docker image is built
1. Once complete, run `docker run -i -p 8888:8888 -t -v {path of local folder where scripts are}:/code/ {name of container}_app /bin/bash` An example is `docker run -i -t  -p 8888:8888 -v /mnt/diatom-pipeline:/code/ diatom-pipeline_app /bin/bash`
2. Now to copy PEAR to the user/bin of the container `cp pear/pear /usr/local/bin`
3. You can now access the notebook by entering `http://localhost:8888` in your browser. 
4. Next, run `sh diatomPipeline.dms sequences lookuptable.txt`
4. Enter `exit`



## Results
After running to completion the pipeline will output two files. The first file, **Abundances.fail.csv**, contains in the top row a list of the IDs of the samples which failed QC,i.e. samples returning <3000 sequences after quality trimming and read merging.

The second file, **Abundances.pass.csv**, contains a list of all samples that have passed QC.Column 1 lists the identity (strain ID) of each of strains identified in each sample and row 1 lists the sample ID (taken from the lookuptable.txtfile. Each cell will show the number of reads mapping to a particular strain for a particular sample.


### Tips and tricks

- You can find the container ID/ name by running `docker ps`.
- Large datasets need lots of RAM to run. The more your machine has, the better. 8 GB is minimum.
- Do not be concerned if the pipeline takes some time to run, our tests have shown that for a year's worth of data over 20 hours is normal.
