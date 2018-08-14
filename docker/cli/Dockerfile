FROM centos:6
MAINTAINER http://github.com/radiomics

RUN yum update -y && \
  yum groupinstall -y "Development Tools" && \
  yum install -y curl \
   curl-devel \
   coreutils \
   freetype-devel \
   gcc \
   gcc-c++ \
   gettext \
   libpng-devel \
   openssl-devel \
   perl \
   perl-CPAN \
   perl-ExtUtils-MakeMaker \
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

# Build and install git from source.
WORKDIR /usr/src
ENV GIT_VERSION 2.5.0
RUN wget --no-check-certificate https://www.kernel.org/pub/software/scm/git/git-${GIT_VERSION}.tar.gz && \
  tar xvzf git-${GIT_VERSION}.tar.gz && \
  cd git-${GIT_VERSION} && \
  ./configure --prefix=/usr && \
  make && \
  make install && \
  cd .. && rm -rf git-${GIT_VERSION}*

# Build and install CMake from source.
WORKDIR /usr/src
RUN git clone git://cmake.org/cmake.git CMake && \
  cd CMake && \
  git checkout v3.7.2 && \
  mkdir /usr/src/CMake-build && \
  cd /usr/src/CMake-build && \
  /usr/src/CMake/bootstrap \
    --parallel=$(grep -c processor /proc/cpuinfo) \
    --prefix=/usr && \
  make -j$(grep -c processor /proc/cpuinfo) && \
  ./bin/cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_USE_OPENSSL:BOOL=ON . && \
  make install && \
  cd .. && rm -rf CMake*

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

# Build and install ninja from source.
WORKDIR /usr/src
RUN git clone https://github.com/martine/ninja.git && \
  cd ninja && \
  git checkout v1.6.0 && \
  ./configure.py --bootstrap && \
  mv ninja /usr/bin/ && \
  cd .. && rm -rf ninja

# Setup necessary python packages
WORKDIR /usr/src
RUN python --version && \
  wget --no-check-certificate https://bootstrap.pypa.io/get-pip.py && \
  python get-pip.py && \
  pip install wheel>=0.29.0 && \
  pip install setuptools>=28.0.0 && \
  pip install scipy && \
  pip install trimesh && \
  pip install scikit-image && \
  pip install --trusted-host www.itk.org -f https://itk.org/SimpleITKDoxygen/html/PyDownloadPage.html SimpleITK>=0.9.1 && \
  python -c "import SimpleITK; print('SimpleITK Version:' + SimpleITK.Version_VersionString())"

# Install PyRadiomics
WORKDIR /usr/src
RUN git clone https://github.com/radiomics/pyradiomics.git && \
  cd pyradiomics && \
  pip install -r requirements.txt && \
  python setup.py install

WORKDIR /usr/src
ENTRYPOINT ["pyradiomics"]
