# Build Pyradiomics inside the Jupyter Datascience Notebook

FROM jupyter/datascience-notebook

MAINTAINER Radiomics Project (http://github.com/radiomics)

# Build information
ARG BUILD_DATE
ARG GIT_REF

LABEL org.label-schema.build-data=$BUILD_DATE \
      org.label-schema.vcs-url="https://github.com/radiomics/pyradiomics.git" \
      org.label-schema.vcs-ref=$GIT_REF \
      org.label-schema.schema-version="1.0.0-rc1"      

USER root
ADD . /root/pyradiomics
# Install in Python 3
RUN /bin/bash -c "source activate root \
    && cd /root/pyradiomics \
    && python -m pip install --no-cache-dir -r requirements.txt \
    && python setup.py install"
# Python 2
RUN /bin/bash -c "source activate python2 \
    && cd /root/pyradiomics \
    && python -m pip install --no-cache-dir -r requirements.txt \
    && python setup.py install"

# Install sample data and notebooks
ADD data/ /home/jovyan/work/example_data/
ADD bin/Notebooks/RadiomicsExample.ipynb /home/jovyan/work/
ADD bin/Notebooks/FeatureVisualization.ipynb /home/jovyan/work/
ADD bin/Notebooks/FeatureVisualizationWithClustering.ipynb /home/jovyan/work/
ADD bin/Notebooks/FilteringEffects.ipynb /home/jovyan/work/
ADD bin/Params.yaml /home/jovyan/work/

# Make a global directory and link it to the work directory
RUN mkdir /data
RUN ln -s /data /home/jovyan/work/data

RUN chown -R jovyan:users /home/jovyan/work

# Trust the notebooks that we've installed
USER jovyan
RUN jupyter trust /home/jovyan/work/*.ipynb

# Run the notebooks
RUN jupyter nbconvert --ExecutePreprocessor.kernel_name=python2 --ExecutePreprocessor.timeout=-1 --to notebook --execute /home/jovyan/work/*.ipynb
RUN jupyter nbconvert --ExecutePreprocessor.kernel_name=python3 --ExecutePreprocessor.timeout=-1 --to notebook --execute /home/jovyan/work/*.ipynb

# The user's data will show up as /data
VOLUME /data
