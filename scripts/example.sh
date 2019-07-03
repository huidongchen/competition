# build the method container
CONTAINER_FOLDER=containers/methods/r/
TAG=dynverse/r_example

chmod +x $CONTAINER_FOLDER/main.R
docker build -t $TAG $CONTAINER_FOLDER

# create a results folder to store the results
RESULT_FOLDER=$(pwd)/results
RESULT_MOUNT="-v $RESULT_FOLDER:/ti/"
mkdir $RESULT_FOLDER

# define the folder where the data is stored
DATA_FOLDER=$(pwd)/examples/inputs
GT_FOLDER=$(pwd)/examples/ground-truths
DATA_MOUNT="-v $DATA_FOLDER:/data/"
GT_MOUNT="-v $GT_FOLDER:/gt"
DATASET_ID=linear

# run the R container on the test dataset
docker run $DATA_MOUNT $RESULT_MOUNT dynverse/r_example /data/${DATASET_ID}.h5 /ti/
ls $RESULT_FOLDER

# convert the output csv's into an HDF5 file that can be used for further processing (https://dynverse.org)
docker run $DATA_MOUNT $RESULT_MOUNT dynverse/convert_output --dataset /data/${DATASET_ID}.h5 --output_folder /ti/ --model /ti/model.h5
ls $RESULT_FOLDER

# then evaluate it using the dyneval docker
docker run $DATA_MOUNT $GT_MOUNT $RESULT_MOUNT dynverse/dyneval --groundtruth /gt/${DATASET_ID}.h5 --model /ti/model.h5 --output_scores /ti/scores.json
cat $RESULT_FOLDER/scores.json
