#!/usr/bin/env python

import argparse
import re
import math
import subprocess
from snakemake.utils import read_job_properties

# Parse command line arguments
parser = argparse.ArgumentParser(description="Submits a Snakemake job to M3")
parser.add_argument("job_script", help="Job script created by Snakemake")
args = parser.parse_args()

job_properties = read_job_properties(args.job_script)

# Create an m3 submission command and submit the job to the cluster
cluster_params = job_properties["params"]["cluster"]
memory = cluster_params["memory"]
time = cluster_params["time"]
cores = cluster_params["cores"]
qos = cluster_params["qos"]
partition = cluster_params["partition"]

#############################################
# EDIT BELOW HERE FOR CLUSTER CUSTOMIZATION #
#############################################
# If multiple cores are requested, then we need to request memory per CPU, not total memory
if cores > 1:
    memory_int = int(re.match("([0-9]*)[^0-9]", memory).group(1))
    memory_units = re.match("[0-9]*([^0-9]*)$", memory).group(1)
    memory_per_cpu = math.floor(float(memory_int) / float(cores))
    memory = str(memory_per_cpu) + memory_units

submission_command = ["sbatch"]
submission_command += ["--job-name=MetaLAFFA", "--ntasks=1", "--qos=%s" % qos, "--partition=%s" % partition, "--mem-per-cpu=%s" % memory, "--time=%s" % time, "--cpus-per-task=%s" % cores]
submission_command += [args.job_script]
subprocess.run(submission_command)
