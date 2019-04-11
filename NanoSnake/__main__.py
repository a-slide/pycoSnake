#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
import os
import sys
import argparse
import pkg_resources
import logging
import shutil
from collections import *
import inspect

# Third party library
from snakemake import snakemake

# Local imports
from NanoSnake import __version__ as package_version
from NanoSnake import __name__ as package_name
from NanoSnake import __description__ as package_description
from NanoSnake.common import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~LOGGING INFO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GLOBAL DIRS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ENVS_DIR = pkg_resources.resource_filename (package_name, "envs")
SCRIPTS_DIR = pkg_resources.resource_filename (package_name, "scripts")
RULES_DIR = pkg_resources.resource_filename (package_name, "rules")
SNAKEFILES_DIR = pkg_resources.resource_filename (package_name, "snakefiles")
TEMPLATES_DIR = pkg_resources.resource_filename (package_name, "templates")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CLI ENTRY POINT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def main(args=None):
    """
    Main entry point for NanoSnake command line interface
    """

    # Parser and subparsers for command
    parser = argparse.ArgumentParser (description=package_description)
    parser.add_argument("--version", action="version", version="{} v{}".format(package_name, package_version))
    subparsers = parser.add_subparsers (description="NanoSnake implements the following subcommands", dest="subcommands")
    subparsers.required = True

    # Template subparser
    subparser_t = subparsers.add_parser("generate_template", description="Generate teample files corresponding to the indicated workflow")
    subparser_t.set_defaults(func=generate_template)
    subparser_t.add_argument('subcommands_name', type=str, choices=['DNA_methylation', 'RNA_counts'], help="Workflow to generate template files for.")
    subparser_t.add_argument("--outdir", default="./", type=str, help="Path where to write the template files (default: %(default)s).")
    subparser_t.add_argument("--overwrite", action="store_true", default=False, help="Overwrite existing files if they already exist (default: %(default)s).")

    # DNA_methylation subparser
    subparser_dm = subparsers.add_parser("DNA_methylation", description="Workflow to evaluate DNA methylation in Nanopore data using Nanopolish")
    subparser_dm.set_defaults(func=DNA_methylation)
    subparser_dm_IO = subparser_dm.add_argument_group("input/output options")
    subparser_dm_IO.add_argument("--reference", "-r", required=True, type=str, help="Path to a Fasta reference file to be used for read mapping (Required if neither `config_template` nor `sample_template` is given)")

    # RNA_counts subparser
    subparser_rc = subparsers.add_parser("RNA_counts", description="Workflow aligning RNA on the transcriptome an generating estimated counts")
    subparser_rc.set_defaults(func=RNA_counts)

    # Add common group parsers
    for sp in [subparser_dm, subparser_rc]:
        sp_IO = add_argument_group (sp, "input/output options")
        sp_IO.add_argument("--sample_sheet", "-s", required=True, type=str, help="Path to a tabulated sample sheet")
        sp_IO.add_argument("--config_file", "-c", default=None, type=str, help="Snakemake configuration YAML file (default: %(default)s)")
        sp_IO.add_argument("--workdir", "-d", default="./", type=str, help="Path to the working dir where to deploy the workflow (default: %(default)s)")
        sp_verbosity = sp.add_mutually_exclusive_group()
        sp_verbosity.add_argument("--verbose", "-v", action="store_true", default=False, help="Show additional debug output (default: %(default)s)")
        sp_verbosity.add_argument("--quiet", "-q", action="store_true", default=False, help="Reduce overall output (default: %(default)s)")
        sp_snakemake = add_argument_group(sp, "Snakemake options")
        sp_snakemake.add_argument("--report", type=str, default=None, help="create an HTML report for a previous run at the given path")
        sp_snakemake.add_argument("--listrules", action="store_true", default=False, help="list rules")
        sp_snakemake.add_argument("--list_target_rules", action="store_true", default=False, help="list target rules")
        sp_snakemake.add_argument("--cores", type=int, default=1, help="the number of provided cores (ignored when using cluster support)")
        sp_snakemake.add_argument("--nodes", type=int, default=1, help="the number of provided cluster nodes (ignored without cluster support)")
        sp_snakemake.add_argument("--local_cores", type=int, default=1, help="the number of provided local cores if in cluster mode (ignored without cluster support)")
        sp_snakemake.add_argument("--targets", type=str, nargs='+', default=[], help="list of targets, e.g. rule or file names")
        sp_snakemake.add_argument("--dryrun", action="store_true", default=False, help="only dry-run the workflow")
        sp_snakemake.add_argument("--touch", action="store_true", default=False, help="only touch all output files if present")
        sp_snakemake.add_argument("--forcetargets", action="store_true", default=False, help="force given targets to be re-created")
        sp_snakemake.add_argument("--forceall", action="store_true", default=False, help="force all output files to be re-created")
        sp_snakemake.add_argument("--forcerun", type=str, nargs='+', default=[], help="list of files and rules that shall be re-created/re-executed")
        sp_snakemake.add_argument("--prioritytargets", type=str, nargs='+', default=[], help="list of targets that shall be run with maximum priority")
        sp_snakemake.add_argument("--stats", type=str, default=None, help="path to file that shall contain stats about the workflow execution")
        sp_snakemake.add_argument("--printreason", action="store_true", default=False, help="print the reason for the execution of each job")
        sp_snakemake.add_argument("--printshellcmds", action="store_true", default=False, help="print the shell command of each job")
        sp_snakemake.add_argument("--printdag", action="store_true", default=False, help="print the dag in the graphviz dot language")
        sp_snakemake.add_argument("--printrulegraph", action="store_true", default=False, help="print the graph of rules in the graphviz dot language")
        sp_snakemake.add_argument("--printd3dag", action="store_true", default=False, help="print a D3.js compatible JSON representation of the DAG")
        sp_snakemake.add_argument("--nocolor", action="store_true", default=False, help="do not print colored output")
        sp_snakemake.add_argument("--keepgoing", action="store_true", default=False, help="keep goind upon errors")
        sp_snakemake.add_argument("--cluster", type=str, default=None, help="submission command of a cluster or batch system to use, e.g. qsub")
        sp_snakemake.add_argument("--cluster_config", type=str, nargs='+', default=[], help="configuration file for cluster options, or list thereof")
        sp_snakemake.add_argument("--cluster_sync", type=str, default=None, help="blocking cluster submission command (like SGE ‘qsub -sync y’)")
        sp_snakemake.add_argument("--drmaa", type=str, default=None, help="if not None use DRMAA for cluster support, str specifies native args passed to the cluster when submitting a job")
        sp_snakemake.add_argument("--drmaa_log_dir", type=str, default=None, help="the path to stdout and stderr output of DRMAA jobs")
        sp_snakemake.add_argument("--jobname", type=str, default='snakejob.{rulename}.{jobid}.sh', help="naming scheme for cluster job scripts")
        sp_snakemake.add_argument("--immediate_submit", action="store_true", default=False, help="immediately submit all cluster jobs, regardless of dependencies")
        sp_snakemake.add_argument("--ignore_ambiguity", action="store_true", default=False, help="ignore ambiguous rules and always take the first possible one")
        sp_snakemake.add_argument("--unlock", action="store_true", default=False, help="just unlock the working directory")
        sp_snakemake.add_argument("--cleanup_metadata", type=str, nargs='+', default=[], help="just cleanup metadata of given list of output files")
        sp_snakemake.add_argument("--cleanup_conda", action="store_true", default=False, help="just cleanup unused conda environments")
        sp_snakemake.add_argument("--cleanup_shadow", action="store_true", default=False, help="just cleanup old shadow directories")
        sp_snakemake.add_argument("--force_incomplete", action="store_true", default=False, help="force the re-creation of incomplete files")
        sp_snakemake.add_argument("--ignore_incomplete", action="store_true", default=False, help="ignore incomplete files")
        sp_snakemake.add_argument("--list_version_changes", action="store_true", default=False, help="list output files with changed rule version")
        sp_snakemake.add_argument("--list_code_changes", action="store_true", default=False, help="list output files with changed rule code")
        sp_snakemake.add_argument("--list_input_changes", action="store_true", default=False, help="list output files with changed input files")
        sp_snakemake.add_argument("--list_params_changes", action="store_true", default=False, help="list output files with changed params")
        sp_snakemake.add_argument("--list_untracked", action="store_true", default=False, help="list files in the workdir that are not used in the workflow")
        sp_snakemake.add_argument("--archive", type=str, default=None, help="archive workflow into the given tarball")
        sp_snakemake.add_argument("--delete_all_output", action="store_true", default=False, help="remove all files generated by the workflow")
        sp_snakemake.add_argument("--delete_temp_output", action="store_true", default=False, help="remove all temporary files generated by the workflow")
        sp_snakemake.add_argument("--latency_wait", type=int, default=3, help="how many seconds to wait for an output file to appear after the execution of a job, e.g. to handle filesystem latency")
        sp_snakemake.add_argument("--wait_for_files", type=str, nargs='+', default=[], help="wait for given files to be present before executing the workflow")
        sp_snakemake.add_argument("--list_resources", action="store_true", default=False, help="list resources used in the workflow")
        sp_snakemake.add_argument("--summary", action="store_true", default=False, help="list summary of all output files and their status. If no option is specified a basic summary will be ouput. If ‘detailed’ is added as an option e.g –summary detailed, extra info about the input and shell commands will be included")
        sp_snakemake.add_argument("--detailed_summary", action="store_true", default=False, help="list summary of all input and output files and their status")
        sp_snakemake.add_argument("--print_compilation", action="store_true", default=False, help="print the compilation of the snakefile")
        sp_snakemake.add_argument("--debug", action="store_true", default=False, help="allow to use the debugger within rules")
        sp_snakemake.add_argument("--notemp", action="store_true", default=False, help="ignore temp file flags, e.g. do not delete output files marked as temp after use")
        sp_snakemake.add_argument("--keep_remote_local", action="store_true", default=False, help="keep local copies of remote files")
        sp_snakemake.add_argument("--nodeps", action="store_true", default=False, help="ignore dependencies")
        sp_snakemake.add_argument("--keep_target_files", action="store_true", default=False, help="do not adjust the paths of given target files relative to the working directory.")
        sp_snakemake.add_argument("--jobscript", type=str, default=None, help="path to a custom shell script template for cluster jobs")
        sp_snakemake.add_argument("--overwrite_shellcmd", type=str, default=None, help="a shell command that shall be executed instead of those given in the workflow. This is for debugging purposes only")
        sp_snakemake.add_argument("--updated_files", type=str, nargs='+', default=[], help="a list that will be filled with the files that are updated or created during the workflow execution")
        sp_snakemake.add_argument("--max_jobs_per_second", type=int, default=None, help="maximal number of cluster/drmaa jobs per second, None to impose no limit")
        sp_snakemake.add_argument("--restart_times", type=int, default=0, help="number of times to restart failing jobs")
        sp_snakemake.add_argument("--attempt", type=int, default=1, help="initial value of Job.attempt. This is intended for internal use only")
        sp_snakemake.add_argument("--force_use_threads", action="store_true", default=False, help="whether to force use of threads over processes. helpful if shared memory is full or unavailable")
        sp_snakemake.add_argument("--use_singularity", action="store_true", default=False, help="run jobs in singularity containers (if defined with singularity directive)")
        sp_snakemake.add_argument("--singularity_args", type=str, default=None, help="additional arguments to pass to singularity")
        sp_snakemake.add_argument("--conda_prefix", type=str, default=None, help="the directory in which conda environments will be created")
        sp_snakemake.add_argument("--singularity_prefix", type=str, default=None, help="the directory to which singularity images will be pulled")
        sp_snakemake.add_argument("--shadow_prefix", type=str, default=None, help="prefix for shadow directories. The job-specific shadow directories will be created in $SHADOW_PREFIX/shadow/")
        sp_snakemake.add_argument("--create_envs_only", action="store_true", default=False, help="if specified, only builds the conda environments specified for each job, then exits.")
        sp_snakemake.add_argument("--list_conda_envs", action="store_true", default=False, help="list conda environments and their location on disk.")
        sp_snakemake.add_argument("--wrapper_prefix", type=str, default=None, help="prefix for wrapper script URLs")
        sp_snakemake.add_argument("--kubernetes", type=str, default=None, help="submit jobs to kubernetes, using the given namespace.")
        sp_snakemake.add_argument("--kubernetes_envvars", type=str, nargs='+', default=[], help="environment variables that shall be passed to kubernetes jobs.")
        sp_snakemake.add_argument("--container_image", type=str, default=None, help="Docker image to use, e.g., for kubernetes.")
        sp_snakemake.add_argument("--default_remote_provider", type=str, default=None, help="default remote provider to use instead of local files (e.g. S3, GS)")
        sp_snakemake.add_argument("--default_remote_prefix", type=str, default=None, help="prefix for default remote provider (e.g. name of the bucket).")
        sp_snakemake.add_argument("--cluster_status", type=str, default=None, help="status command for cluster execution. If None, Snakemake will rely on flag files. Otherwise, it expects the command to return “success”, “failure” or “running” when executing with a cluster jobid as single argument.")
        sp_snakemake.add_argument("--export_cwl", type=str, default=None, help="Compile workflow to CWL and save to given file")
        #sp_snakemake.add_argument("--use_conda", action="store_true", default=False, help="create conda environments for each job (defined with conda directive of rules)")

    # Parse args and call subfunction
    args = parser.parse_args()

    # Alter logging level if required
    if args.subcommands != "generate_template":
        if args.verbose:
            logger.setLevel (logging.DEBUG)
        elif args.quiet:
            logger.setLevel (logging.WARNING)

    args.func(args)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SUBPARSERS FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def generate_template (args):
    """"""
    for fn in ("samples.tsv", "config.yaml"):
        dest_file = os.path.join(args.outdir, fn)
        if os.path.isfile(dest_file) and not args.overwrite:
            raise NanoSnakeError (f"Template file {dest_file} already exists. Please use --overwrite if you want to replace the existing file")
        else:
            src_file = os.path.join (TEMPLATES_DIR, "{}_{}".format(args.subcommands_name, fn))
            logger.info (f"Create template file {dest_file} ")
            shutil.copyfile (src_file, dest_file)

def DNA_methylation (args):
    """"""
    # Define default package files and dir
    snake_fn = os.path.join (SNAKEFILES_DIR, "DNA_methylation_snakefile.py")

    # Add additional config options
    logger.info ("Build config dict for snakemake")
    config = {
        "sample_sheet": args.sample_sheet,
        "reference": args.reference,
        "envs_dir": ENVS_DIR,
        "scripts_dir": SCRIPTS_DIR,
        "rules_dir": RULES_DIR}
    logger.debug (config)

    # Default config if not given
    if args.config_file:
        config_file = args.config_file
    else:
        config_file = os.path.join (TEMPLATES_DIR, "DNA_methylation_config.yaml")

    # Filter other args option compatible with snakemake API
    kwargs = filter_valid_snakemake_options (args)
    logger.debug (kwargs)

    # Run Snakemake through the API
    logger.warning ("RUNING SNAKEMAKE PIPELINE")
    snakemake (snakefile=snake_fn, configfile=config_file, config=config, use_conda=True, **kwargs)

def RNA_counts (args):
    """"""
    pass

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HELPER FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def add_argument_group (parser, title):
    """Add group only is it doesn't exist yet"""
    for group in parser._action_groups:
        if group.title == title:
            return group
    return parser.add_argument_group(title)

def filter_valid_snakemake_options (args):
    """Filter out options that are not in the snakemake API"""
    valid_options = list(inspect.signature(snakemake).parameters.keys())
    valid_kwargs = OrderedDict()
    for k,v in vars(args).items():
        if k in valid_options:
            valid_kwargs[k] = v
    return valid_kwargs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SCRIPT ENTRY POINT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == "__main__":
    main ()
