FROM radiomics/pyradiomics:CLI
MAINTAINER https://github.com/Radiomics

# Install required dcm2niix converions tools
WORKDIR /usr/src
RUN git clone https://github.com/rordenlab/dcm2niix.git && \
  cd dcm2niix && \
  mkdir build && cd build && \
  cmake .. && \
  make install

# Install additional needed python packages
# Prevent pandas from installing version 0.25.0, as this version is unable to import
# module _lzma (in package lmza). See also https://github.com/pandas-dev/pandas/issues/27543
WORKDIR /usr/src
RUN pip3 install pandas!=0.25.0 && \
  pip3 install pydicom

# Copy the required DCMQI tools from the latest DCMQI docker image
COPY --from=qiicr/dcmqi:latest /usr/bin/tid1500writer /usr/bin/segimage2itkimage /usr/bin/
# Copy the pyradiomics-dcm.py script to ensure latest changes are included in the image
COPY ./labs/pyradiomics-dcm/* /usr/src/pyradiomics/labs/pyradiomics-dcm/

WORKDIR /usr/src/pyradiomics/labs/pyradiomics-dcm
ENTRYPOINT ["python3", "pyradiomics-dcm.py"]
