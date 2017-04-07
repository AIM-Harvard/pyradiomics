# Build Pyradiomics inside the Jupyter Datascience Notebook

FROM jupyter/datascience-notebook

MAINTAINER Radiomics Project (http://github.com/radiomics)

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
ADD bin/Notebooks/FilteringEffects.ipynb /home/jovyan/work/
ADD bin/Params.yaml /home/jovyan/work/

# Make a global directory and link it to the work directory
RUN mkdir /data
RUN ln -s /data /home/jovyan/work/data

RUN chown -R jovyan:users /home/jovyan/work

# Trust the notebooks that we've installed
USER jovyan
RUN jupyter trust /home/jovyan/work/RadiomicsExample.ipynb
RUN jupyter trust /home/jovyan/work/FeatureVisualization.ipynb
RUN jupyter trust /home/jovyan/work/FilteringEffects.ipynb

# The user's data will show up as /data
VOLUME /data
