vpkg_require anaconda/2.5.0:python2
vpkg_require openmpi/1.10.2-intel64-2016
rm -f PythonOutput.txt
python Debug.py > PythonOutput.txt
