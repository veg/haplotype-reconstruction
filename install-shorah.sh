#  This script will install ShoRAH into your current directory.
#  to run it, type 
#
#    bash install-shorah.sh
#
#  Virtualenv needs to be installed prior to running this script.
#  A new directory named 'shorah' will be created.
#  ShoRAH's python script library will be located in ./shorah/bin
#  after running this script.

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

configCommand="./configure --prefix="$shorahDir" PYTHON="$shorahDir"/bin/python"

eval "$configCommand"

make -j4

make install

printf "\n\n\n          ShoRAH 1.1.3 has been installed.\n\n\n\n"
