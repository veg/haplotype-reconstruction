#          
#        VPhaser/VProfiler install script
# 
#	 This script will install VPhaser and VProfiler in your current directory.
#  To run it, type:
#
#    bash install-vphaser.sh
#
#  The zipped package file from the Broad Institute is required to run this script.
#  See the accompanying file vphaser-info.txt for more information.
#


conda config --add channels bioconda

printf "\nChecking dependencies. This may take a few minutes.\n\n"

printf "Checking dependencies (1/6)...\n"
conda install -c defaults -y perl

printf "Checking dependencies (2/6)...\n"
conda install -c r -y r

printf "Checking dependencies (3/6)...\n"
conda install -c r -y r-gplots

printf "Checking dependencies (4/6)...\n"
conda install -c bioconda -y samtools

printf "Checking dependencies (5/6)...\n"
conda install -c bioconda -y mosaik

printf "Checking dependencies (6/6)...\n"
conda install -c bioconda -y muscle

unzip v_phaser.zip

mv __MACOSX VpSoftwarePackage

mv VpSoftwarePackage vphaser-vprofiler

cd vphaser-vprofiler

printf "Configuring VPhaser and VProfiler...\n"

perl configPaths.pl configfile.txt

printf "\nInstallation complete. Running sample data to verify installation...\n\n"

mkdir installtest

cp -R TestData/* installtest/

cd installtest

printf "Validating VPhaser..."

perl ../vphaser.pl -i rcVTest_final.qlx -o vp_VTest

printf "Validating VProfiler...\n"

perl ../vprofiler.pl -i vprofiler_input_VTest.txt -o vpro -noendvariant=10 -nt -codon -haplo -haploseq

if diff -q vpro_haplotypes/ ExpectedResults/vpro_haplotypes/; then
	printf "\n\n      VPhaser and VProfiler installed successfully.\n\n\n"
else
	printf "Warning: test output does not match sample output.\n"
fi

cd ..

rm -r installtest
