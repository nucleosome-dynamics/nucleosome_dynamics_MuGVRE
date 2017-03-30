#!/usr/bin/env sh

bin="/home/rilla/nucleServ/nucleosomeDynamics_wf.py"
root_dir="/home/rilla/nucleServ/wrap_test"
public_dir="/orozco/services/Rdata/MuG/MuG_public/"
config="/home/rilla/nucleServ/wrap_test/config.json"
metadata="/home/rilla/nucleServ/wrap_test/metadata.json"
out_metadata="/home/rilla/nucleServ/wrap_test/out_metadata.json"

python3 $bin                 \
    --config $config         \
    --root_dir $root_dir     \
    --public_dir $public_dir \
    --metadata $metadata     \
    --out_metadata $out_metadata
