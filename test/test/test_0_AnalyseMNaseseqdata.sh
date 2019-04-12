#!/bin/bash

###
### Testing in a local installation
### the VRE server CMD
###
### * Automatically created by MuGVRE *
###


# Local installation - EDIT IF REQUIRED

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

TOOL_EXECUTABLE=$CWD/../../nucleosome_dynamics.py # Main tool wrapper
INPUTS_DIR=$CWD/../data                           # Input data directory 
TEST_DATA_DIR=$CWD/json/0_AnalyseMNaseseqdata     # Test configuration files (JSONs)
WORKING_DIR=$CWD/../example_run                   # Define output directory



# Adapting test configuration files to INPUTS_DIR

sed -e 's#$WORKING_DIR#'$WORKING_DIR'#' $TEST_DATA_DIR/config.json.template      > $TEST_DATA_DIR/config.json
sed -e 's#$INPUTS_DIR#'$INPUTS_DIR'#'   $TEST_DATA_DIR/in_metadata.json.template > $TEST_DATA_DIR/in_metadata.json


# Running ND2 tool ...

if [ -d  $WORKING_DIR ]; then rm -r $WORKING_DIR/; mkdir -p $WORKING_DIR; else mkdir -p $WORKING_DIR; fi
cd $WORKING_DIR

echo "--- Test execution: $WORKING_DIR"
echo "--- Start time: `date`"

time $TOOL_EXECUTABLE --config $TEST_DATA_DIR/config.json --in_metadata $TEST_DATA_DIR/in_metadata.json --out_metadata $WORKING_DIR/out_metadata.json > $WORKING_DIR/tool.log

