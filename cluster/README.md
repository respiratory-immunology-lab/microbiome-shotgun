## Working on the cluster

Here we provide basic commands for connecting and submitting jobs on the cluster. See https://docs.massive.org.au/ for more information on the cluster.

```
# Connect to cluster
ssh [username]@m3.massive.org.au # Mac users
ssh -l [username] m3.massive.org.au # Linux users

# Start a new smux session (n) for example with 20 cpus (--ntasks) for 1 day (--time)
smux n --ntasks=20 --time=1-00:00:00

# List ongoing smux sessions
smux l

# List ongoing jobs
show_job

# Load modules
module load [module]

# Cancel a smux session
scancel [jobID]
```

## Downloading data from BaseSpace

See https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview for downloading data from BaseSpace straight to the cluster. 

```
# Download data from BaseSpace
./bs -c Australia download project -i [projectID] -o [directory]
```
