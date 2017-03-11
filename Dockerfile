FROM centos:5
MAINTAINER http://github.com/radiomics

RUN yum update -y && \
  yum groupinstall -y "Development Tools" && \
  yum install -y curl \
   curl-devel \
   coreutils \
   gcc \
   gcc-c++ \
   gettext \
   openssl-devel \
   perl \
   wget \
   zlib-devel

WORKDIR /etc/yum.repos.d
RUN wget http://people.centos.org/tru/devtools-2/devtools-2.repo
RUN yum install -y devtoolset-2-gcc \
  devtoolset-2-binutils \
  devtoolset-2-gcc-gfortran \
  devtoolset-2-gcc-c++
ENV CC /opt/rh/devtoolset-2/root/usr/bin/gcc
ENV CXX /opt/rh/devtoolset-2/root/usr/bin/g++
ENV FC /opt/rh/devtoolset-2/root/usr/bin/gfortran

# Build and install git from source.
#WORKDIR /usr/src
#ENV GIT_VERSION 2.11.0
#RUN wget https://github.com/git/git/archive/v${GIT_VERSION}.tar.gz -O git-${GIT_VERSION}.tar.gz && \
#  tar xvzf git-${GIT_VERSION}.tar.gz && \
#  cd git-${GIT_VERSION} && \
#  make prefix=/usr all && \
#  make prefix=/usr install && \
#  cd .. && rm -rf git-${GIT_VERSION}*

# Build and install Python from source.
WORKDIR /usr/src
ENV PYTHON_VERSION 2.7.10
RUN wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz && \
  tar xvzf Python-${PYTHON_VERSION}.tgz && \
  cd Python-${PYTHON_VERSION} && \
  ./configure && \
  make -j$(grep -c processor /proc/cpuinfo) && \
  make install && \
  cd .. && rm -rf Python-${PYTHON_VERSION}*

# Install pyradiomics
WORKDIR /usr/src
ADD . pyradiomics/
RUN cd pyradiomics && \
  python --version && \
  wget --no-check-certificate https://bootstrap.pypa.io/get-pip.py && \
  python get-pip.py && \
  pip install wheel>=0.29.0 && \
  pip install setuptools>=28.0.0 && \
  pip install --trusted-host www.itk.org -f https://itk.org/SimpleITKDoxygen/html/PyDownloadPage.html SimpleITK>=0.9.1 && \
  python -c "import SimpleITK; print('SimpleITK Version:' + SimpleITK.Version_VersionString())" && \
  pip install -r requirements.txt && \
  python setup.py install

WORKDIR /usr/src
ENTRYPOINT ["pyradiomics"]
