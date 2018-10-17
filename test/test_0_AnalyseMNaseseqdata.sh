#!/bin/bash

###
### Testing in a local installation
### the VRE server CMD
###
### * Automatically created by MuGVRE *
###

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# Local installation - EDIT IF REQUIRED

WORKING_DIR=$CWD/test_run001
TOOL_EXECUTABLE=$CWD/../nucleosome_dynamics.py


# Input files for test run - ADAPT TO YOUR OWN DATA IF REQUIRED

TEST_JSON_DIR=$CWD/json/0_AnalyseMNaseseqdata
TEST_DATA_DIR=$CWD/data/

sed -i -e 's#"name": "execution", "value":.*$#"name": "execution", "value": "'"$WORKING_DIR"'"#' $TEST_JSON_DIR/config.json
sed -i -e 's#"file_path": .*cellcycleG2_chrII.bam"#"file_path": "'"$TEST_DATA_DIR"'cellcycleG2_chrII.bam"#' $TEST_JSON_DIR/in_metadata.json
sed -i -e 's#"file_path": .*cellcycleG2_chrII.bam.bai"#"file_path": "'"$TEST_DATA_DIR"'cellcycleG2_chrII.bam.bai"#' $TEST_JSON_DIR/in_metadata.json
sed -i -e 's#"file_path": .*cellcycleM_chrII.bam"#"file_path": "'"$TEST_DATA_DIR"'cellcycleM_chrII.bam"#' $TEST_JSON_DIR/in_metadata.json
sed -i -e 's#"file_path": .*cellcycleM_chrII.bam.bai"#"file_path": "'"$TEST_DATA_DIR"'cellcycleM_chrII.bam.bai"#' $TEST_JSON_DIR/in_metadata.json
sed -i -e 's#"file_path": .*refGenomes.*,#"file_path": "'"$TEST_DATA_DIR"'refGenomes/",#' $TEST_JSON_DIR/in_metadata.json



# Running Nucleosome Dynamics full workflow ...

if [ -d  $WORKING_DIR ]; then rm -r $WORKING_DIR/; mkdir -p $WORKING_DIR; else mkdir -p $WORKING_DIR; fi
cd $WORKING_DIR

COMMAND="$TOOL_EXECUTABLE --config $TEST_JSON_DIR/config.json --in_metadata $TEST_JSON_DIR/in_metadata.json --out_metadata $WORKING_DIR/out_metadata.json > $WORKING_DIR/tool.log"

echo "--- Test execution: $WORKING_DIR"
echo "--- Start time: `date`"
echo "--- $COMMAND"
echo ""

eval time $COMMAND
