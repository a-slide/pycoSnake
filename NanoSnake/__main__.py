#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
import os
import sys
import argparse
import pkg_resources
from collections import *
import glob

# Third party library
from snakemake import snakemake
from loguru import logger as log

# Local imports
from NanoSnake import __version__ as package_version
from NanoSnake import __name__ as package_name
from NanoSnake import __description__ as package_description
from NanoSnake.common import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GLOBAL DIRS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
WORKFLOW_DIR = pkg_resources.resource_filename (package_name, "workflows")
DATA_DIR = pkg_resources.resource_filename (package_name, "test_data")
WRAPPER_DIR = pkg_resources.resource_filename (package_name, "wrappers")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~LIST EXISTING WRAPPERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
WRAPPER_PREFIX = "file:{}/".format(WRAPPER_DIR)
WRAPPERS = []
for w in glob.glob(os.path.join(WRAPPER_DIR, "*")):
    if os.path.isfile(os.path.join(w, "snakefile.py")):
        WRAPPERS.append(os.path.split(w)[-1])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CLI ENTRY POINT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def main(args=None):
    """ Main entry point for NanoSnake command line interface """

    # Parser and subparsers for command
    parser = argparse.ArgumentParser (description=package_description)
    parser.add_argument("--version", action="version", version="{} v{}".format(package_name, package_version))
    subparsers = parser.add_subparsers (description="%(prog)s implements the following subcommands", dest="subcommands")
    subparsers.required = True

    # test_wrappers subparser
    subparser_tw = subparsers.add_parser("test_wrappers", description="Test Nanosnake wrappers")
    subparser_tw.set_defaults(func=test_wrappers, type="test")
    subparser_tw.add_argument("--wrappers", "-w", default=WRAPPERS, nargs='+', choices=WRAPPERS, type=str, help="List of wrappers to test (default: all)")
    subparser_tw.add_argument("--keep_output", "-k", action="store_true", default=False, help="Keep temporary output files generated during tests (default: %(default)s)")
    subparser_tw.add_argument("--clean_output", "-c", action="store_true", default=False, help="clean all temporary output files generated during tests (default: %(default)s)")
    subparser_tw.add_argument("--cores", "-j", type=int, default=1, help="the number of provided cores (default: %(default)s)")
    subparser_tw.add_argument("--workdir", "-d", default="./", type=str, help="Path to the working dir where to deploy the workflow (default: %(default)s)")

    # DNA_ONT subparser
    subparser_dna_ont = subparsers.add_parser("DNA_ONT", description="Workflow for DNA Analysis of Nanopore data")
    subparser_dna_ont.set_defaults(func=DNA_ONT, type="workflow")
    subparser_dna_ont_IO = subparser_dna_ont.add_argument_group("input/output options")
    subparser_dna_ont_IO.add_argument("--genome", "-g", default=None, type=str, help="Path to an ENSEMBL FASTA reference genome file/URL to be used for read mapping (required)")
    subparser_dna_ont_IO.add_argument("--annotation", "-a", default=None, type=str, help="Path to an ENSEMBL GFF3 annotation file/URL containing transcript annotations (required)")
    subparser_dna_ont_IO.add_argument("--sample_sheet", "-s", default=None, type=str, help="Path to a tabulated sample sheet (required)")

    # RNA_illumina subparser
    subparser_rna_illumina = subparsers.add_parser("RNA_illumina", description="Workflow for RNA Analysis of Illumina data")
    subparser_rna_illumina.set_defaults(func=RNA_illumina, type="workflow")
    subparser_rna_illumina_IO = subparser_rna_illumina.add_argument_group("input/output options")
    subparser_rna_illumina_IO.add_argument("--genome", "-g", default=None, type=str, help="Path to an ENSEMBL FASTA reference genome file/URL to be used for read mapping (required)")
    subparser_rna_illumina_IO.add_argument("--transcriptome", "-t", default=None, type=str, help="Path to an ENSEMBL cDNA FASTA reference transcriptome file/URL to be used for read counting (required)")
    subparser_rna_illumina_IO.add_argument("--annotation", "-a", default=None, type=str, help="Path to an ENSEMBL GFF3 annotation file/URL containing transcript annotations (required)")
    subparser_rna_illumina_IO.add_argument("--sample_sheet", "-s", default=None, type=str, help="Path to a tabulated sample sheet (required)")

    # pycoMeth subparser
    subparser_pycometh = subparsers.add_parser("pycoMeth", description="Workflow for DNA Analysis of Nanopore data")
    subparser_pycometh.set_defaults(func=pycoMeth, type="workflow")
    subparser_pycometh_IO = subparser_pycometh.add_argument_group("input/output options")
    subparser_pycometh_IO.add_argument("--genome", "-g", default=None, type=str, help="Path to an ENSEMBL FASTA reference genome file/URL to be used for read mapping (required)")
    subparser_pycometh_IO.add_argument("--annotation", "-a", default=None, type=str, help="Path to an ENSEMBL GFF3 annotation file/URL containing transcript annotations (required)")
    subparser_pycometh_IO.add_argument("--sample_sheet", "-s", default=None, type=str, help="Path to a tabulated sample sheet (required)")
    subparser_pycometh_INT = subparser_pycometh.add_argument_group("Interval aggregation options")
    subparser_pycometh_INT.add_argument("--interval_mode", "-i", nargs="?", default="cpg_islands", choices=["cpg_islands", "external_bed", "sliding_window"],
        help="""How to aggregate CpG into intervals. (default: %(default)s)
            * cpg_islands: Find CpG islands in the genome file.
            * external_bed: use intervals provided in an external bed file.
            * sliding_window: Use a sliding window along entire genome.""")
    subparser_pycometh_INT.add_argument("--external_bed", "-b", default=None, type=str, help="Path to a bed file containing intervals, if interval_mode is set to external_bed")

    # Add common options for all parsers
    for sp in [subparser_dna_ont, subparser_rna_illumina, subparser_pycometh, subparser_tw]:
        sp_verbosity = sp.add_mutually_exclusive_group()
        sp_verbosity.add_argument("--verbose", "-v", action="store_true", default=False, help="Show additional debug output (default: %(default)s)")
        sp_verbosity.add_argument("--quiet", "-q", action="store_true", default=False, help="Reduce overall output (default: %(default)s)")

    # Add common options for workflow parsers
    for sp in [subparser_dna_ont, subparser_rna_illumina, subparser_pycometh]:
        sp_IO = add_argument_group (sp, "input/output options")
        sp_IO.add_argument("--config", "-c", default=None, type=str, help="Snakemake configuration YAML file (required in local mode)")
        sp_IO.add_argument("--cluster_config", default=None, type=str, help="Snakemake cluster configuration YAML file (required in cluster mode)")
        sp_IO.add_argument("--workdir", "-d", default="./", type=str, help="Path to the working dir where to deploy the workflow (default: %(default)s)")
        sp_template = add_argument_group (sp, "Template options")
        sp_template.add_argument("--generate_template", "-e", type=str, nargs="+", default=[], choices=["all", "sample_sheet", "config", "cluster_config"], help="Generate template files (configs + sample_sheet) in workdir and exit (default: %(default)s)")
        sp_template.add_argument("--overwrite_template", "-o", action="store_true", default=False, help="Overwrite existing template files if they already exist (default: %(default)s)")
        sp_snakemake = add_argument_group(sp, "Snakemake options")
        sp_snakemake.add_argument("--report", type=str, default=None, help="create an HTML report for a previous run at the given path (default: %(default)s)")
        sp_snakemake.add_argument("--listrules", action="store_true", default=False, help="list rules (default: %(default)s)")
        sp_snakemake.add_argument("--list_target_rules", action="store_true", default=False, help="list target rules (default: %(default)s)")
        sp_snakemake.add_argument("--cores", "-j", type=int, default=1, help="the number of provided cores (default: %(default)s)")
        sp_snakemake.add_argument("--nodes", type=int, default=1, help="the number of provided cluster nodes (ignored without cluster support) (default: %(default)s)")
        sp_snakemake.add_argument("--targets", type=str, nargs='+', default=[], help="list of targets, e.g. rule or file names (default: %(default)s)")
        sp_snakemake.add_argument("--dryrun", action="store_true", default=False, help="only dry-run the workflow (default: %(default)s)")
        sp_snakemake.add_argument("--touch", action="store_true", default=False, help="only touch all output files if present (default: %(default)s)")
        sp_snakemake.add_argument("--forcetargets", action="store_true", default=False, help="force given targets to be re-created (default: %(default)s)")
        sp_snakemake.add_argument("--forceall", action="store_true", default=False, help="force all output files to be re-created (default: %(default)s)")
        sp_snakemake.add_argument("--forcerun", type=str, nargs='+', default=[], help="list of files and rules that shall be re-created/re-executed (default: %(default)s)")
        sp_snakemake.add_argument("--prioritytargets", type=str, nargs='+', default=[], help="list of targets that shall be run with maximum priority (default: %(default)s)")
        sp_snakemake.add_argument("--stats", type=str, default=None, help="path to file that shall contain stats about the workflow execution (default: %(default)s)")
        sp_snakemake.add_argument("--printreason", action="store_true", default=False, help="print the reason for the execution of each job (default: %(default)s)")
        sp_snakemake.add_argument("--printshellcmds", action="store_true", default=False, help="print the shell command of each job (default: %(default)s)")
        sp_snakemake.add_argument("--printdag", action="store_true", default=False, help="print the dag in the graphviz dot language (default: %(default)s)")
        sp_snakemake.add_argument("--printrulegraph", action="store_true", default=False, help="print the graph of rules in the graphviz dot language (default: %(default)s)")
        sp_snakemake.add_argument("--printd3dag", action="store_true", default=False, help="print a D3.js compatible JSON representation of the DAG (default: %(default)s)")
        sp_snakemake.add_argument("--nocolor", action="store_true", default=False, help="do not print colored output (default: %(default)s)")
        sp_snakemake.add_argument("--keepgoing", action="store_true", default=False, help="keep going upon errors (default: %(default)s)")
        sp_snakemake.add_argument("--cluster", type=str, default=None, help="submission command of a cluster or batch system to use, e.g. qsub (default: %(default)s)")
        sp_snakemake.add_argument("--drmaa_log_dir", type=str, default=None, help="the path to stdout and stderr output of DRMAA jobs (default: %(default)s)")
        sp_snakemake.add_argument("--jobname", type=str, default='snakejob.{rulename}.{jobid}.sh', help="naming scheme for cluster job scripts (default: %(default)s)")
        sp_snakemake.add_argument("--immediate_submit", action="store_true", default=False, help="immediately submit all cluster jobs, regardless of dependencies (default: %(default)s)")
        sp_snakemake.add_argument("--ignore_ambiguity", action="store_true", default=False, help="ignore ambiguous rules and always take the first possible one (default: %(default)s)")
        sp_snakemake.add_argument("--unlock", action="store_true", default=False, help="just unlock the working directory (default: %(default)s)")
        sp_snakemake.add_argument("--cleanup_metadata", type=str, nargs='+', default=[], help="just cleanup metadata of given list of output files (default: %(default)s)")
        sp_snakemake.add_argument("--cleanup_conda", action="store_true", default=False, help="just cleanup unused conda environments (default: %(default)s)")
        sp_snakemake.add_argument("--force_incomplete", action="store_true", default=False, help="force the re-creation of incomplete files (default: %(default)s)")
        sp_snakemake.add_argument("--ignore_incomplete", action="store_true", default=False, help="ignore incomplete files (default: %(default)s)")
        sp_snakemake.add_argument("--archive", type=str, default=None, help="archive workflow into the given tarball (default: %(default)s)")
        sp_snakemake.add_argument("--delete_all_output", action="store_true", default=False, help="remove all files generated by the workflow (default: %(default)s)")
        sp_snakemake.add_argument("--delete_temp_output", action="store_true", default=False, help="remove all temporary files generated by the workflow (default: %(default)s)")
        sp_snakemake.add_argument("--latency_wait", type=int, default=3, help="how many seconds to wait for an output file to appear after the execution of a job, e.g. to handle filesystem latency (default: %(default)s)")
        sp_snakemake.add_argument("--wait_for_files", type=str, nargs='+', default=[], help="wait for given files to be present before executing the workflow (default: %(default)s)")
        sp_snakemake.add_argument("--list_resources", action="store_true", default=False, help="list resources used in the workflow (default: %(default)s)")
        sp_snakemake.add_argument("--summary", action="store_true", default=False, help="list summary of all output files and their status. If no option is specified a basic summary will be ouput. If ‘detailed’ is added as an option e.g –summary detailed, extra info about the input and shell commands will be included (default: %(default)s)")
        sp_snakemake.add_argument("--detailed_summary", action="store_true", default=False, help="list summary of all input and output files and their status (default: %(default)s)")
        sp_snakemake.add_argument("--print_compilation", action="store_true", default=False, help="print the compilation of the snakefile (default: %(default)s)")
        sp_snakemake.add_argument("--debug", action="store_true", default=False, help="allow to use the debugger within rules (default: %(default)s)")
        sp_snakemake.add_argument("--notemp", action="store_true", default=False, help="ignore temp file flags, e.g. do not delete output files marked as temp after use (default: %(default)s)")
        sp_snakemake.add_argument("--keep_remote_local", action="store_true", default=False, help="keep local copies of remote files (default: %(default)s)")
        sp_snakemake.add_argument("--nodeps", action="store_true", default=False, help="ignore dependencies (default: %(default)s)")
        sp_snakemake.add_argument("--keep_target_files", action="store_true", default=False, help="do not adjust the paths of given target files relative to the working directory. (default: %(default)s)")
        sp_snakemake.add_argument("--overwrite_shellcmd", type=str, default=None, help="a shell command that shall be executed instead of those given in the workflow. This is for debugging purposes only (default: %(default)s)")
        sp_snakemake.add_argument("--updated_files", type=str, nargs='+', default=[], help="a list that will be filled with the files that are updated or created during the workflow execution (default: %(default)s)")
        sp_snakemake.add_argument("--max_jobs_per_second", type=int, default=None, help="maximal number of cluster/drmaa jobs per second, None to impose no limit (default: %(default)s)")
        sp_snakemake.add_argument("--restart_times", type=int, default=0, help="number of times to restart failing jobs (default: %(default)s)")
        sp_snakemake.add_argument("--attempt", type=int, default=1, help="initial value of Job.attempt. This is intended for internal use only (default: %(default)s)")
        sp_snakemake.add_argument("--force_use_threads", action="store_true", default=False, help="whether to force use of threads over processes. helpful if shared memory is full or unavailable (default: %(default)s)")
        sp_snakemake.add_argument("--conda_prefix", type=str, default=None, help="the directory in which conda environments will be created (default: %(default)s)")
        sp_snakemake.add_argument("--create_envs_only", action="store_true", default=False, help="if specified, only builds the conda environments specified for each job, then exits. (default: %(default)s)")
        sp_snakemake.add_argument("--list_conda_envs", action="store_true", default=False, help="list conda environments and their location on disk. (default: %(default)s)")
        sp_snakemake.add_argument("--wrapper_prefix",type=str, default=WRAPPER_PREFIX, help=" (prefix for wrapper script URLs. (default: %(default)s)")
#        sp_snakemake.add_argument("--use_singularity", action="store_true", default=False, help="run jobs in singularity containers (if defined with singularity directive) (default: %(default)s)")
#        sp_snakemake.add_argument("--singularity_args", type=str, default=None, help="additional arguments to pass to singularity (default: %(default)s)")
#        sp_snakemake.add_argument("--singularity_prefix", type=str, default=None, help="the directory to which singularity images will be pulled (default: %(default)s)")
#        sp_snakemake.add_argument("--drmaa", type=str, default=None, help="if not None use DRMAA for cluster support, str specifies native args passed to the cluster when submitting a job (default: %(default)s)")
#        sp_snakemake.add_argument("--cluster_sync", type=str, default=None, help="blocking cluster submission command (like SGE ‘qsub -sync y’) (default: %(default)s)")
#        sp_snakemake.add_argument("--cleanup_shadow", action="store_true", default=False, help="just cleanup old shadow directories (default: %(default)s)")
#        sp_snakemake.add_argument("--list_version_changes", action="store_true", default=False, help="list output files with changed rule version (default: %(default)s)")
#        sp_snakemake.add_argument("--list_code_changes", action="store_true", default=False, help="list output files with changed rule code (default: %(default)s)")
#        sp_snakemake.add_argument("--list_input_changes", action="store_true", default=False, help="list output files with changed input files (default: %(default)s)")
#        sp_snakemake.add_argument("--list_params_changes", action="store_true", default=False, help="list output files with changed params (default: %(default)s)")
#        sp_snakemake.add_argument("--list_untracked", action="store_true", default=False, help="list files in the workdir that are not used in the workflow (default: %(default)s)")
#        sp_snakemake.add_argument("--jobscript", type=str, default=None, help="path to a custom shell script template for cluster jobs (default: %(default)s)")
#        sp_snakemake.add_argument("--shadow_prefix", type=str, default=None, help="prefix for shadow directories. The job-specific shadow directories will be created in $SHADOW_PREFIX/shadow/ (default: %(default)s)")
#        sp_snakemake.add_argument("--kubernetes", type=str, default=None, help="submit jobs to kubernetes, using the given namespace. (default: %(default)s)")
#        sp_snakemake.add_argument("--kubernetes_envvars", type=str, nargs='+', default=[], help="environment variables that shall be passed to kubernetes jobs. (default: %(default)s)")
#        sp_snakemake.add_argument("--container_image", type=str, default=None, help="Docker image to use, e.g., for kubernetes. (default: %(default)s)")
#        sp_snakemake.add_argument("--default_remote_provider", type=str, default=None, help="default remote provider to use instead of local files (e.g. S3, GS) (default: %(default)s)")
#        sp_snakemake.add_argument("--default_remote_prefix", type=str, default=None, help="prefix for default remote provider (e.g. name of the bucket). (default: %(default)s)")
#        sp_snakemake.add_argument("--cluster_status", type=str, default=None, help="status command for cluster execution. If None, Snakemake will rely on flag files. Otherwise, it expects the command to return “success”, “failure” or “running” when executing with a cluster jobid as single argument. (default: %(default)s)")
#        sp_snakemake.add_argument("--export_cwl", type=str, default=None, help="Compile workflow to CWL and save to given file (default: %(default)s)")

    # Parse args and and define logger verbose level
    args = parser.parse_args()
    set_log_level(quiet=args.quiet, verbose=args.verbose)
    log.warning ("RUNNING {} v{}".format(package_name, package_version))

    if args.type == "workflow":

        # Unlock locked dir and exit
        if args.unlock:
            log.warning (f"Unlocking working directory")
            unlock_dir (workdir=args.workdir)
            sys.exit()

        # Generate templates and exit
        if args.generate_template:
            log.warning (f"Generate template files in working directory")
            generate_template (
                workflow_dir=WORKFLOW_DIR,
                templates=args.generate_template,
                workflow=args.subcommands,
                workdir=args.workdir,
                overwrite=args.overwrite_template,
                verbose=args.verbose,
                quiet=args.quiet)
            sys.exit()

        # Cluster stuff to simplify options
        if args.cluster_config:
            log.warning (f"INITIALISING WORKFLOW IN CLUSTER MODE")
            args.local_cores = get_yaml_val(yaml_fn=args.cluster_config, val_name="cluster_cores", default=args.cores)
            args.nodes = get_yaml_val(yaml_fn=args.cluster_config, val_name="cluster_nodes", default=args.nodes)
            args.cluster = get_yaml_val(yaml_fn=args.cluster_config, val_name="cluster_cmd", default=args.cluster)
            args.config = args.cluster_config
            log.debug (f"Cores:{args.local_cores} / Nodes:{args.nodes} / Cluster_cmd:{args.cluster}")
        elif args.config :
            log.warning (f"INITIALISING WORKFLOW IN LOCAL MODE")
        else:
            raise NanoSnakeError("A configuration file `--config` or a cluster configuration file `--cluster_config` is required")

    # Run workflow
    args.func(args)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DNA_ONT SUBPARSER FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def DNA_ONT (args):
    """"""
    # Get and check config files
    log.warning ("CHECKING CONFIGURATION FILES")
    snakefile = get_snakefile_fn(workflow_dir=WORKFLOW_DIR, workflow=args.subcommands)
    configfile = get_config_fn(config=args.config)

    # Store additionnal options to pass to snakemake
    log.info ("Build config dict for snakemake")
    config = {
        "genome":required_option("genome", args.genome),
        "annotation":required_option("annotation", args.annotation),
        "sample_sheet":get_sample_sheet(sample_sheet=args.sample_sheet, required_fields=["sample_id", "fastq", "fast5", "seq_summary"])}
    log.debug (config)

    # Filter other args option compatible with snakemake API
    kwargs = filter_valid_snakemake_options (args)
    log.debug (kwargs)

    # Run Snakemake through the API
    log.warning ("RUNNING SNAKEMAKE PIPELINE")
    snakemake (snakefile=snakefile, configfile=configfile, config=config, use_conda=True, **kwargs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RNA_illumina SUBPARSER FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def RNA_illumina (args):
    """"""
    # Get and check config files
    log.warning ("CHECKING CONFIGURATION FILES")
    snakefile = get_snakefile_fn(workflow_dir=WORKFLOW_DIR, workflow=args.subcommands)
    configfile = get_config_fn(config=args.config)

    # Store additionnal options to pass to snakemake
    log.info ("Build config dict for snakemake")
    config = {
        "genome":required_option("genome", args.genome),
        "transcriptome":required_option("transcriptome", args.transcriptome),
        "annotation":required_option("annotation", args.annotation),
        "sample_sheet":get_sample_sheet(sample_sheet=args.sample_sheet, required_fields=["sample_id", "fastq1", "fastq2"])}
    log.debug (config)

    # Filter other args option compatible with snakemake API
    kwargs = filter_valid_snakemake_options (args)
    log.debug (kwargs)

    # Run Snakemake through the API
    log.warning ("RUNNING SNAKEMAKE PIPELINE")
    snakemake (snakefile=snakefile, configfile=configfile, config=config, use_conda=True, **kwargs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TEST SUBPARSER FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def pycoMeth (args):
    """"""
    # Get and check config files
    log.warning ("CHECKING CONFIGURATION FILES")
    snakefile = get_snakefile_fn(workflow_dir=WORKFLOW_DIR, workflow=args.subcommands)
    configfile = get_config_fn(config=args.config)

    # Store additionnal options to pass to snakemake
    log.info ("Build config dict for snakemake")
    config = {
        "genome":required_option("genome", args.genome),
        "annotation":required_option("annotation", args.annotation),
        "sample_sheet":get_sample_sheet(sample_sheet=args.sample_sheet, required_fields=["sample_id", "methylation_calls"]),
        "interval_mode":args.interval_mode,
        "external_bed":args.external_bed}
    log.debug (config)

    # Filter other args option compatible with snakemake API
    kwargs = filter_valid_snakemake_options (args)
    log.debug (kwargs)

    # Run Snakemake through the API
    log.warning ("RUNNING SNAKEMAKE PIPELINE")
    snakemake (snakefile=snakefile, configfile=configfile, config=config, use_conda=True, **kwargs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TEST SUBPARSER FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def test_wrappers (args):
    """"""
    # Cleanup data and leave
    if args.clean_output:
        log.info("Removing output data")
        for wrapper_name in args.wrappers:
            wrapper_workdir = os.path.join(args.workdir, wrapper_name)
            shutil.rmtree(wrapper_workdir, ignore_errors=True)
        sys.exit()

    # Test wrappers
    for wrapper_name in args.wrappers:
        log.warning("Testing Wrapper {}".format(wrapper_name))
        try:
            snakefile = get_snakefile_fn(workflow_dir=WRAPPER_DIR, workflow=wrapper_name)
            wrapper_workdir = os.path.join(args.workdir, wrapper_name)
            log.debug("Working in directory: {}".format(wrapper_workdir))

            #Run Snakemake through the API
            snakemake (
                snakefile = snakefile,
                workdir = wrapper_workdir,
                config = {"data_dir":DATA_DIR},
                wrapper_prefix = WRAPPER_PREFIX,
                use_conda = True,
                cores = args.cores,
                verbose = args.verbose,
                quiet = args.quiet)

        finally:
            log.debug("List of file generated: {}".format(os.listdir(wrapper_workdir)))
            shutil.rmtree(os.path.join(wrapper_workdir, ".snakemake"), ignore_errors=True)
            if not args.keep_output:
                log.debug("Removing temporary directory")
                shutil.rmtree(wrapper_workdir, ignore_errors=True)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Helper functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def filter_valid_snakemake_options (args):
    """Filter out options that are not in the snakemake API"""
    valid_options = list(inspect.signature(snakemake).parameters.keys())
    valid_kwargs = OrderedDict()
    for k,v in vars(args).items():
        if k in valid_options and k != "config":
            valid_kwargs[k] = v
    return valid_kwargs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SCRIPT ENTRY POINT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == "__main__":
    main ()
