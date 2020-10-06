# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:latest

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3-pip \
    gnupg \
    mongodb



# install miniconda
ENV MINICONDA_VERSION 4.8.2
ENV CONDA_DIR ~/miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_$MINICONDA_VERSION-Linux-x86_64.sh -O ~/miniconda.sh && \
  chmod +x ~/miniconda.sh && \
  ~/miniconda.sh -b -p $CONDA_DIR && \
  rm ~/miniconda.sh


# make non-activate conda commands available
ENV PATH=/root/miniconda3/bin:$PATH

# make conda activate command available from /bin/bash --login shells
RUN echo ". /root/miniconda3/etc/profile.d/conda.sh" >> ~/.profile

# RUN conda info --envs

# build the conda environment
RUN mkdir /envs
ENV ENV_PREFIX /envs/mongo_env
RUN conda update --name base --channel defaults conda && \
    conda create -p $ENV_PREFIX -c conda-forge mongo-tools --force && \
    conda clean --all --yes

# Make RUN commands use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]

# make conda activate command available from /bin/bash --interative shells
RUN conda init bash

#RUN wget https://www.mongodb.org/static/pgp/server-4.4.asc
#RUN sudo apt-key add server-4.4.asc

#RUN sh -c 'wget -qO - https://www.mongodb.org/static/pgp/server-4.4.asc | sudo #apt-key add -'

RUN echo "deb [ arch=amd64,arm64 ] https://repo.mongodb.org/apt/ubuntu focal/mongodb-org/4.4 multiverse" >> /etc/apt/sources.list.d/mongodb-org-4.4.list


## update system libraries
#RUN apt-get update
#&& \
    #apt-get upgrade -y && \
    #apt-get clean

#RUN apt-get update && \
        #apt-get install -y --force-yes pwgen mongodb-org-server

#RUN sudo apt-get install -y mongodb-org

# copy necessary files
## app folder
COPY app.R ./app/
COPY requirements.txt ./app/
COPY data ./data
## renv.lock file
#COPY /example-app/renv.lock ./renv.lock

RUN conda activate /envs/mongo_env

RUN mkdir /data/db

# RUN mongod --fork --dbpath /data/db --logpath /app/mongod.log
WORKDIR /data
RUN mongod --fork --dbpath /data/db --logpath /app/mongod.log && \
    /envs/mongo_env/bin/mongoimport -d asp_bgc_fa -c info --file info_fa.json && \
    /envs/mongo_env/bin/mongoimport -d asp_bgc_fa -c bgc --file bgc_fa.json && \
    mongod --shutdown &&\
    echo "OK"

RUN echo "mongod --fork --dbpath /data/db --logpath /app/mongod.log" >> ~/.profile

#RUN pip3 install --upgrade pip
#RUN pip install --no-cache-dir --force-reinstall -r app/requirements.txt

# install renv & restore packages
# RUN Rscript -e 'install.packages("renv")'
# RUN Rscript -e 'renv::restore()'

RUN Rscript -e  'install.packages("tidyverse")'
RUN Rscript -e  'install.packages("gggenes")'
RUN Rscript -e  'install.packages("RColorBrewer")'
RUN Rscript -e  'install.packages("RCurl")'
RUN Rscript -e  'install.packages("DT")'
RUN Rscript -e  'install.packages("UpSetR")'
RUN Rscript -e  'install.packages("mongolite")'
RUN Rscript -e  'install.packages("shinythemes")'
RUN Rscript -e  'install.packages("shinycssloaders")'
RUN Rscript -e  'install.packages("gggenomes")'




# expose port
EXPOSE 3838

# run app on container start
CMD ["mongod", "--fork", "--dbpath", "/data/db", "--logpath", "/app/mongod3.log", "&&", "R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]
