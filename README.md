Transform ascii models into a single A3d files. More capabilities will be added later

### Step 0 - Load or create the conda environment
conda create env -f environment.yml
conda activate <conda_env_name>

### Step 1 - Create the sample files
Modify conf.yaml to adjust to your own paths
./0_create_sample.sh <path_to_ascii_files>
This will output 3 files that you can name in conf.yaml, useful for the next step

### Step 2 - Run the conversion script
python 1_ascii2A3d.py conf.yaml
This step will output the a3d file. Right now it does only one component at a time, which will require you to stitch a3d files together (copy and paste) if you have multiple parameters (let's say Vs and Xi) but this will be updated in the future.
