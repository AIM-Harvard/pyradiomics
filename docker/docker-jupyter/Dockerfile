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

# Jupyter datascience notebook does not support python 2 anymore, install it manually.
# "source activate python2 || (...)" first checks if the python 2 environment exists,
# and only installs if the environment does not exist
# See also https://github.com/jupyter/docker-stacks/issues/432
# Next, also install python 2 kernel globally, so it can be found from the root
RUN /bin/bash -c "source activate python2 \
    || (conda create -n python2 python=2 ipykernel \
    && pip install kernda --no-cache \
    && $CONDA_DIR/envs/python2/bin/python -m ipykernel install \
    && kernda -o -y /usr/local/share/jupyter/kernels/python2/kernel.json \
    && pip uninstall kernda -y)"

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
ADD data/ /home/jovyan/work/data/
ADD notebooks/RadiomicsExample.ipynb /home/jovyan/work/notebooks/
ADD notebooks/FeatureVisualization.ipynb /home/jovyan/work/notebooks/
ADD notebooks/FeatureVisualizationWithClustering.ipynb /home/jovyan/work/notebooks/
ADD notebooks/FilteringEffects.ipynb /home/jovyan/work/notebooks/
ADD notebooks/helloRadiomics.ipynb /home/jovyan/work/notebooks/
ADD notebooks/helloFeatureClass.ipynb /home/jovyan/work/notebooks/
ADD notebooks/PyRadiomicsExample.ipynb /home/jovyan/work/notebooks/
ADD examples/exampleSettings/Params.yaml /home/jovyan/work/examples/exampleSettings/

# Make a global directory and link it to the work directory
RUN mkdir /data
RUN ln -s /data /home/jovyan/work/data

RUN chown -R jovyan:users /home/jovyan/work

# Trust the notebooks that we've installed
USER jovyan
RUN jupyter trust /home/jovyan/work/notebooks/*.ipynb

# Run the notebooks
RUN jupyter nbconvert --ExecutePreprocessor.kernel_name=python2 --ExecutePreprocessor.timeout=-1 --to notebook --execute /home/jovyan/work/notebooks/helloRadiomics.ipynb /home/jovyan/work/notebooks/helloFeatureClass.ipynb /home/jovyan/work/notebooks/PyRadiomicsExample.ipynb
RUN jupyter nbconvert --ExecutePreprocessor.kernel_name=python3 --ExecutePreprocessor.timeout=-1 --to notebook --execute /home/jovyan/work/notebooks/helloRadiomics.ipynb /home/jovyan/work/notebooks/helloFeatureClass.ipynb /home/jovyan/work/notebooks/PyRadiomicsExample.ipynb

# The user's data will show up as /data
VOLUME /data
