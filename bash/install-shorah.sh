#  This script will install ShoRAH into your current directory.
#  to run it, type
#
#    bash install-shorah.sh
#
#  The program virtualenv needs to be installed prior to running this script.
#  A new directory named 'shorah' will be created.
#  ShoRAH's python script library will be created in ./shorah/bin
#
#  NOTE: If you want to install ShoRAH using conda, you can do so by typing:
#
#    conda install -c bioconda -y shorah

virtualenv shorah

cd shorah

installDir=$(pwd)

source bin/activate

pip install numpy Biopython

wget https://github.com/cbg-ethz/shorah/releases/download/v1.1.3/shorah-1.1.3.tar.bz2

bunzip2 shorah-1.1.3.tar.bz2

tar -xvf shorah-1.1.3.tar

rm shorah-1.1.3.tar

cd shorah-1.1.3

configCommand="./configure --prefix="$installDir" PYTHON="$installDir"/bin/python"

eval "$configCommand"

make -j4

make install

printf "\n  ShoRAH 1.1.3 install script finished.\n\n"
