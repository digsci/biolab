FROM ubuntu

MAINTAINER Lancaster University version 1.0 <t.byrne1@lancaster.ac.uk>

# get build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils && \
    echo "installed ubuntu" && \
    apt install -y zlibc zlib1g zlib1g-dev &&  \
    echo "installed zlibc" && \
    apt-get install -y python-qt4 && \
    apt install -y curl && \
    apt install -y  moreutils && \
    apt install -y parallel && \
    apt-get install -y  build-essential autoconf automake libtool &&  \
    apt-get install -y libbz2-dev && \
    apt-get install -y  wget &&  \
    apt install -y libgl1-mesa-glx && \
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
    sh Miniconda2-latest-Linux-x86_64.sh -p /miniconda -b
    ENV PATH=/miniconda/bin:${PATH}
    RUN conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda && \
    conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --add channels r && \
    conda config --add channels bioconda && \
    conda install -y biopython &&  \
    conda install -y -c bioconda multiqc multiqc=1.5 && \
    conda install -y sickle-trim && \
    conda install -c bioconda blast=2.2 && \
    apt-get install -y libxml2-dev libxslt1-dev && \
    apt install -y  python-dev && \
    pip install lxml==3.5.0 &&  \
    pip install pygal &&  \
    pip install numpy==1.10.4 && \
    pip install pandas==0.19.0 && \
    pip install cutadapt==1.9.1 && \
    pip install qiime==1.9.1

	
    #now to install jupyter

    RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py" && \
    python3 get-pip.py --user

    RUN pip install jupyter

    # This will be our application root folder
    RUN mkdir /code
    WORKDIR /code


    CMD jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root



    # Create a new system user
    RUN useradd -ms /bin/bash jupyter

    ADD ./pear/pear /usr/local/bin

    # Change to this new user
    USER jupyter
