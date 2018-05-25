installDir=$(pwd)

virtualenv shorah

cd shorah

#echo "$installDir"

shorahDir=$(pwd)

source bin/activate

pip install numpy Biopython

wget https://github.com/cbg-ethz/shorah/releases/download/v1.1.3/shorah-1.1.3.tar.bz2

bunzip2 shorah-1.1.3.tar.bz2

tar -xvf shorah-1.1.3.tar

rm shorah-1.1.3.tar

cd shorah-1.1.3

configDir=$(pwd)

configCommand="./configure --prefix="$shorahDir" PYTHON="$shorahDir"/bin/python

make -j4

make install





























