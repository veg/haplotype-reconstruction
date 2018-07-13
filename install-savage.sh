conda create --name savage -y python=2.7

source activate savage

printf "\n\nInstalling dependencies (1/1).."

conda install -c bioconda -y biopython

printf "\n\nInstalling SAVAGE.."

conda install -c bioconda -y savage

mkdir savage-tests
