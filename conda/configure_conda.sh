conda config --set always_yes yes --set changeps1 no
conda config --set anaconda_upload no
conda config --add channels simpleitk
conda config --add channels conda-forge
conda install conda-build
conda install anaconda-client
conda update -q conda
