# Downsample QC

Samples downsampled to 2M reads for QC

Mapping and QC using nf-core/sarek pipeline

## Run commands
Run in screen

Move into separate directory
```
mkdir sarek_GRCh38
cd sarek_GRCh38
```

Setup env
```
source /vulpes/ngi/production/latest/conf/sourceme_sthlm.sh
source activate NGI
```

Run nextflow

```
nextflow run /vulpes/ngi/production/v24.07/sw/sarek/3_4_2/ -profile uppmax --project ngi2016004 -c /vulpes/ngi/production/v24.07/conf/sarek_sthlm.config -c ../nextflow.config -params-file <params.file>
```

Use params.file from below

- human_params.yaml
- mouse_params.yaml
