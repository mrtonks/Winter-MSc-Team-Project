FROM ubuntu


ENV TZ=Europe/Minsk
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get -qq update \
        && apt-get install --no-install-recommends -y \
            libcurl4-openssl-dev \ 
            libssl-dev \
            apt-transport-https \
            python-dev \
            libc-dev \
            python-pip \
            pkg-config \
            liblzma-dev \
            libbz2-dev \
            libpcre3-dev \
            build-essential \
            libblas-dev \
            liblapack-dev \
            gfortran \
            libzmq3-dev \
            curl \
            libfreetype6-dev \
            libpng-dev \
            net-tools \
            procps \
            libreadline-dev \
        && pip install distribute --upgrade \
        && pip install \
            cython \
            pysam \ 
            pyvcf \
        && apt-get autoremove -y \
        && apt-get clean \
	&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN pip install matplotlib

RUN apt-get -qq update \
    && apt-get install -y --no-install-recommends \
        littler \
        r-base \
        r-base-dev \
        r-recommended \
        && echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
        && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
    && ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
    && ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
    && ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
    && ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
    && install.r docopt \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*


RUN install.r \
    doParallel \
    foreach \
    devtools \
    mcclust 
    
RUN installGithub.r \
    keyuan/ccube \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN install.r \
    colorspace \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN mkdir /home/pipeline

COPY ./PyDP-0.2.2.tar.gz /home/pipeline/PyDP-0.2.2.tar.gz
COPY ./pyclone.zip /home/pipeline/pyclone.zip

RUN pip install pyyaml \
    && pip install pandas \
    && cd /home/pipeline/ \
    && tar xvfz /home/pipeline/PyDP-0.2.2.tar.gz \
    && cd /home/pipeline/PyDP-0.2.2/ \
    && python setup.py install \
    && cd /home/pipeline/ \
    && unzip /home/pipeline/pyclone.zip \
    && cd /home/pipeline/pyclone-master/ \
    && python setup.py install \
    && cd /home/

RUN apt-get update && apt-get install -y --no-install-recommends \
	build-essential \
	git \
	gsl-bin \
	libgsl0-dev \
	libnss-sss \
	python-lxml \
	python-numpy \
	python-pip \
	python-qt4 \
	python-scipy \
	python-setuptools \
	python-six

RUN easy_install -U ete2

COPY ./Pipeline /home/pipeline/Pipeline

RUN apt-get install python-tk -y
RUN pip install seaborn
RUN pip install scikit-learn==0.19.0






