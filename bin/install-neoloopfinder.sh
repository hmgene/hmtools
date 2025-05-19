conda config --add channels r
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
mamba create -n neoloop cooler matplotlib pyensembl pybigwig intervaltree rpy2 r-mgcv scikit-learn=1.1.2 joblib=1.1.0 "pomegranate<=0.14.8"
mamba activate neoloop
pip install -U neoloop TADLib
