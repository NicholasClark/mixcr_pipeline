# Running the pipeline

SSH into the HPC and start a tmux session.
```
ssh rhino02
tmux new -s <session_name>
```

Load the Nextflow module
```
module load Nextflow
```
Make a new directory for the run
```
mkdir my_nextflow_run
cd my_nextflow_run
```
Edit the parameters in the tirtl_seq_mixcr_pipeline.nf file:
- Set `params.mixcr_jar_file` to the location of your MiXCR jar file.
- Set `params.output_dir_root` to the folder you want to send output to. The pipeline will make a new subfolder using "--project_name" from your call. 

Run the pipeline:
```
nextflow run <path_to_dir>/tirtl_seq_mixcr_pipeline.nf --input_dir <path_to_fastq_dir> --project_name <your_project_name> -resume
```

The pipeline will run MiXCR on the fastq files for each well (for paired-end sequencing, there are two files per well).
It will create two sub-folders in the output folder: "output_mixcr" and "clones". It also writes a csv file with the list of fastq files found.

All of the MiXCR output for each well, including QC reports, go into the "output_mixcr" folder. The "good" MiXCR output text files with clone counts for each well are moved to the "clones" folder and organized into subfolders by their plate/PCRI-primer (UD01, UD02, etc.). 

Note: PCRI primers are different for alpha/beta (there is a separate "TCRa_UD01" and "TCRb_UD01"), but MiXCR does not know this and will map some of the reads from these primers to the opposite chain. These mismatched files are left in the "output_mixcr" folder and not used.

To save disk space, the pipeline automatically compresses (gzips) text file output from MiXCR.
MiXCR also outputs large ".clns" and ".vdjca" files. We generally delete these to save space.
You may delete these files afterward or un-comment (remove the '#') the lines `rm '${well_id}i'*.vdjca` and `rm '${well_id}i'*.clns` in `tirtl_seq_mixcr_pipeline.nf` to automatically delete these files.

After successfully running MiXCR for all wells, the pipeline runs a custom Python script that
creates an Excel file (.xlsx) with a summary of output for all of the wells and plates/PCRI-primers.

# After the pipeline (TCR-alpha/beta pairing)

After a successful run, you will have a "clones" folder with output folders "UD01", "UD02", etc. for each plate. In each plate's folder, you will have text files (gzipped) of single-chain TCR-alpha and TCR-beta clone counts for each well on the plate.

The rest of the processing, to aggregate the single-chain counts and infer paired TCR-alpha/beta clones for individual samples, which are often one-half or part of a plate, can be done using our "TIRTLtools" R package: https://nicholasclark.github.io/TIRTLtools/.


After installing this package, to get paired TCRs for a sample, you can run the `run_pairing()` function, supplying the directory for the plate (e.g. `your_output_folder/clones/UD01`) and the set of wells used for the sample. Note: wells with low counts for TCR-alpha or TCR-beta will be removed.

See the reference page (https://nicholasclark.github.io/TIRTLtools/reference/run_pairing.html) and the vignette (https://nicholasclark.github.io/TIRTLtools/articles/pairing_updated.html) for more information about pairing.

The pairing function will output three tab-separated (.tsv) files:
- `<sample-name>_TIRTLoutput.tsv(.gz)` – Paired TCRs – a table of computationally paired alpha/beta TCRs
- `<sample-name>_pseudobulk_TRA.tsv(.gz)` – Alpha chain “pseudo-bulk” data – a table of read count summary metrics for alpha-chains across all wells on a plate
- `<sample-name>_pseudobulk_TRB.tsv(.gz)` – Beta chain “pseudo-bulk” data – a table of read count summary metrics for beta-chains across all wells on a plate

## After pairing

See the following vignettes for loading paired data in R  and analyzing it using TIRTLtools:
- https://nicholasclark.github.io/TIRTLtools/articles/data_import.html
- https://nicholasclark.github.io/TIRTLtools/articles/plot_examples.html
