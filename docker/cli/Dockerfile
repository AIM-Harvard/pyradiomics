FROM centos:6
MAINTAINER http://github.com/radiomics

RUN yum update -y && \
    yum groupinstall -y "Development Tools" && \
    yum install -y bzip2-devel \
    curl \
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
RUN git clone https://gitlab.kitware.com/cmake/cmake.git CMake && \
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
ENV PYTHON_VERSION 3.6.9
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

# Setup Python package manager (pip)
WORKDIR /usr/src
RUN python3 --version && \
    wget --no-check-certificate https://bootstrap.pypa.io/get-pip.py && \
    python3 get-pip.py && rm get-pip.py

# Setup necessary python packages
RUN pip3 install wheel>=0.29.0 && \
    pip3 install setuptools>=28.0.0 && \
    pip3 install scipy && \
    pip3 install trimesh && \
    pip3 install scikit-image && \
    pip3 install --trusted-host www.itk.org -f https://itk.org/SimpleITKDoxygen/html/PyDownloadPage.html SimpleITK>=0.9.1 && \
    python3 -c "import SimpleITK; print('SimpleITK Version:' + SimpleITK.Version_VersionString())"

# Install PyRadiomics, then remove the source code to reduce image size
WORKDIR /usr/src
ADD . /usr/src/pyradiomics
RUN cd pyradiomics && \
    pip3 install --no-cache-dir -r requirements.txt && \
    python3 setup.py install && \
    cd .. && rm -rf pyradiomics

WORKDIR /usr/src
ENTRYPOINT ["pyradiomics"]
