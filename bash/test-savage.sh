######################################
#                                    #
#       SAVAGE Testing Script        #
#   For usage, see savage-info.txt   #
#                                    #
######################################

workingDir=$(pwd)

source activate savage

while getopts "hs:" option
do
  case $option in
    h) printf "\n\nUsage: bash test-savage.sh [OPTIONS] [dataset]\n"
			 printf "Example: bash test-savage.sh -s 10 5\n"
			 printf "For more information, see savage-info.txt\n\n"
			 exit 0
    ;;
    s) server=$OPTARG
    ;;
    ?) echo "I didn't recognize that option."
    ;;
  esac
  #shift $((OPTIND-1))
done

dataset=${!OPTIND}

printf "\nDataset: $dataset\n" #debug

if [ -d "savage-tests/$dataset" ]; then
  printf "\nThe directory $workingDir/savage-tests/$dataset already exists. This will overwrite the data in that directory.\n"
  printf "Continue? [y/n]: "
  read response

  if [[ "$response" == "n" || "$response" == "N" ]]; then
    printf "\nExiting SAVAGE test script.\n\n"
    exit
  elif [[ "$response" == "y" || "$response" == "Y" ]]; then
    printf "\nClearing savage-tests/$dataset..."
    eval "rm -r savage-tests/$dataset"
  else
    printf "\nI didn't recognize that response. \
            \nExiting SAVAGE test script.\n\n"
    exit
  fi
fi

eval "mkdir savage-tests/$dataset"

firstDataDir=$(ls -t test-data/454/$dataset/ | head -n 1)

fasta=$(ls -t test-data/454/$dataset/$firstDataDir/BP/1/*Reads.fna | xargs basename)

qual=$(ls -t test-data/454/$dataset/$firstDataDir/BP/1/*Reads.qual | xargs basename)

eval "cp test-data/454/$dataset/$firstDataDir/BP/1/$fasta \
         test-data/454/$dataset/$firstDataDir/BP/1/$qual savage-tests/$dataset"

eval "cd savage-tests/$dataset"

printf "\nCreating FASTQ file...\n"

eval "python $workingDir/fastaqual_to_fastq.py $fasta $qual reads.fastq"


if [ $OPTIND -gt 1 ]; then # Server selected
	printf "\nRunning SAVAGE on server $server. This may take a while.\n"
	eval "time bpsh $server savage --split 1 -s reads.fastq | tee output.txt"
else #No server selected
	printf "\nRunning SAVAGE. This may take a while.\n"
	eval "time savage --split 1 -s reads.fastq | tee output.txt"
fi

testDir=$(pwd)

printf "\nResults have been stored at $testDir/configs_stage_c.fasta\n\n"



