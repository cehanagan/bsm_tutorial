# This will be the Dockerfile for bsm tutorial
# 2 parts: build spotl, build the notebook dependencies

# Build SPOTL 
FROM centos:7 AS compile-spotl
  
RUN yum -y upgrade \
    && yum -y group install "Development Tools" \
    && yum -y install wget vim

WORKDIR /opt

RUN cd /opt/ \
    && wget http://igppweb.ucsd.edu/~agnew/Spotl/spotl.tar.gz \
    && tar -xzf spotl.tar.gz \
    && rm spotl.tar.gz \
    && cd spotl/src/ \
    && echo "FTN = gfortran" >> tempfile \
    && echo "FFLAGS = -O3 -Wuninitialized -fno-f2c -fno-automatic -fno-range-check -fno-backslash" >> tempfile \
    && echo "CC = /usr/bin/gcc" >> tempfile \
    && echo "CFLAGS = -c" >> tempfile \
    && cat Makefile >> tempfile \
    && mv tempfile Makefile \
    && cd ../ \
    && ./install.comp \
    && cd tidmod/ascii \
    && sed -i 's/gzcat/gunzip -c/g' Tobinary \
    && cd ../../ \
    && sed -i 's/Tobinary/.\/Tobinary/g' install.rest \
    && sed -i 's/gzcat/gunzip -c/g' install.rest \
    && sed -i 's/rm /rm -f /g' install.rest \
    && ./install.rest

RUN cd /opt/spotl/working/Exampl/ \
    && sed -i 's/spotl//g' ex1.scr \
    && sed -i 's/spotl//g' ex2.scr \
    && sed -i 's/spotl//g' ex3.scr \
    && sed -i 's/spotl//g' ex4.scr \
    && sed -i 's/spotl//g' ex5.scr \
    && sed -i 's/spotl//g' ex6.scr \
    && sed -i 's/ex5.scr/.\/ex5.scr/g' ex5.scr

# Multi-build the notebook image with the spotl image
# Start from the jupyter/scipy-notebook image, which includes....
# python3, scipy, pandas, matplotlib, ipywidgets (activated), numba, hdf5, etc...
# For more information... https://hub.docker.com/r/jupyter/scipy-notebook/dockerfile

FROM jupyter/scipy-notebook

WORKDIR /home/jovyan/bsm_tutorial

# Copy the build artifact from the SPOTL build stage to this stage
COPY --from=compile-spotl /opt/spotl/ /opt/spotl/

MAINTAINER = UNAVCO Inc. 

USER root

# Change user to default jupyter user (from jupyter/scipy-notebook Dockerfile)
USER $NB_UID

# Install packages with conda
# using channel conda-forge, say yes to all, avoid updating already installed packages
# added sklearn
RUN conda install -c conda-forge --yes --freeze-installed \
        numpy \
        jupyter_contrib_nbextensions \
        obspy \
	xmltodict \
	scikit-learn \
	&& conda clean -afy 

# Install to run SPOTL programs
# (separately so it won't be deleted)
RUN conda install -c conda-forge -y \
	libgfortran

# Copy the notebooks directory
COPY notebooks/* /home/jovyan/bsm_tutorial/notebooks/

# Enable extensions, automatically trust the notebooks
RUN jupyter nbextension enable codefolding/main
RUN jupyter nbextension enable init_cell/main
RUN jupyter nbextension enable jupyter-matplotlib/extension
RUN jupyter trust /home/jovyan/bsm_tutorial/notebooks/NB1.ipynb
RUN jupyter trust /home/jovyan/bsm_tutorial/notebooks/NB2.ipynb
RUN jupyter trust /home/jovyan/bsm_tutorial/notebooks/NB3.ipynb
RUN jupyter trust /home/jovyan/bsm_tutorial/notebooks/NB4.ipynb

# Add SPOTL programs and fortran library paths
ENV PATH="/opt/spotl/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"

# Start the jupyter notebook upon executuion
ENTRYPOINT ["jupyter", "notebook", "--no-browser","--NotebookApp.token=bsm", "--ip=*", "--allow-root", "--notebook-dir=/home/jovyan/bsm_tutorial/notebooks/"]
