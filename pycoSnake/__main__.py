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
import inspect
import textwrap

# Third party library
from snakemake import snakemake
from snakemake import __version__ as snakemake_version
from snakemake.logging import logger, setup_logger

# Local imports
from pycoSnake import __version__ as package_version
from pycoSnake import __name__ as package_name
from pycoSnake import __description__ as package_description
from pycoSnake import workflows_info
from pycoSnake.common import *

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
    """ Main entry point for pycoSnake command line interface """

    # Parser and subparsers for command
    parser = argparse.ArgumentParser (description=package_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--version", action="version", version="{} v{}".format(package_name, package_version))
    description = ["%(prog)s implements the following subcommands"]
    for workflow_name, workflow_info in workflows_info.items():
        description.append("\t* {} v{}: {}".format(workflow_name, workflow_info["version"], workflow_info["description"]))
    subparsers = parser.add_subparsers (description="\n".join(description), dest="subcommand")
    subparsers.required = True

    # test_wrappers subparser
    workflow_name = "test_wrappers"
    workflow_info = workflows_info[workflow_name]
    description = "TEST {} v{}. {}".format(workflow_name, workflow_info["version"], workflow_info["description"])
    subparser_tw = subparsers.add_parser(workflow_name, description=description)
    subparser_tw.set_defaults(parser_func=test_wrappers, workflow_version=workflow_info["version"])
    subparser_tw.add_argument("--wrappers", "-w", default=WRAPPERS, nargs='+', choices=WRAPPERS, type=str, help="List of wrappers to test (default: all)")
    subparser_tw.add_argument("--keep_output", "-k", action="store_true", default=False, help="Keep temporary output files generated during tests (default: %(default)s)")
    subparser_tw.add_argument("--clean_output", "-c", action="store_true", default=False, help="clean all temporary output files generated during tests (default: %(default)s)")
    subparser_tw.add_argument("--cores", "-j", type=int, default=1, help="the number of provided cores (default: %(default)s)")
    subparser_tw.add_argument("--workdir", "-d", default="./", type=str, help="Path to the working dir where to deploy the workflow (default: %(default)s)")

    # DNA_ONT subparser
    workflow_name = "DNA_ONT"
    workflow_info = workflows_info[workflow_name]
    description = [
        "WORKFLOW {} v{}. {}".format(workflow_name, workflow_info["version"], workflow_info["description"]),
        "requires a configuration file and accept any options listed in the official snakemake API documentation:",
        "(https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html)"]
    subparser_dna_ont = subparsers.add_parser(workflow_name, description="\n".join(description))
    subparser_dna_ont.set_defaults(parser_func=workflow_parser, workflow_version=workflow_info["version"])

    # RNA_illumina subparser
    workflow_name = "RNA_illumina"
    workflow_info = workflows_info[workflow_name]
    description = [
        "WORKFLOW {} v{}. {}".format(workflow_name, workflow_info["version"], workflow_info["description"]),
        "requires a configuration file and accept any options listed in the official snakemake API documentation:",
        "(https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html)"]
    subparser_rna_illumina = subparsers.add_parser(workflow_name, description="\n".join(description))
    subparser_rna_illumina.set_defaults(parser_func=workflow_parser, workflow_version=workflow_info["version"])

    # Add common options for workflow parsers
    for sp in [subparser_dna_ont, subparser_rna_illumina]:
        sp.add_argument("--config", "-c", default=None, type=str, help="Snakemake configuration YAML file (required in local mode)")
        sp.add_argument("--cluster_config", default=None, type=str, help="Snakemake cluster configuration YAML file (required in cluster mode)")
        sp.add_argument("--workdir", "-d", default="./", type=str, help="Path to the working dir where to deploy the workflow (default: %(default)s)")
        sp.add_argument("--generate_template", "-e", type=str, nargs="+", default=[], choices=["all", "sample_sheet", "config", "cluster_config"], help="Generate template files (configs + sample_sheet) in workdir and exit (default: %(default)s)")
        sp.add_argument("--overwrite_template", "-o", action="store_true", default=False, help="Overwrite existing template files if they already exist (default: %(default)s)")

    # Add common options for all parsers
    for sp in [subparser_dna_ont, subparser_rna_illumina, subparser_tw]:
        sp_verbosity = sp.add_mutually_exclusive_group()
        sp_verbosity.add_argument("--verbose", "-v", action="store_true", default=False, help="Show additional debug output (default: %(default)s)")
        sp_verbosity.add_argument("--quiet", "-q", action="store_true", default=False, help="Reduce overall output (default: %(default)s)")

    # Parse args and and define logger verbose level
    args, extra = parser.parse_known_args()
    args_dict = autobuild_args(args, extra)

    setup_logger(quiet=args_dict["quiet"], debug=args_dict["verbose"])
    logger.warning ("RUNNING {} v{} with snakemake version v{}".format(package_name, package_version, snakemake_version))
    logger.warning ("RUNNING SUBCOMMAND {} v{}".format(args_dict["subcommand"], args_dict["workflow_version"]))
    args_dict["parser_func"](args_dict)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DNA_ONT SUBPARSER FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def workflow_parser (args_dict):
    """"""
    # Unlock locked dir and exit
    if args_dict.get("unlock", False):
        logger.warning ("Unlocking working directory")
        unlock_dir (workdir=args_dict["workdir"])
        sys.exit()

    # Generate templates and exit
    if args_dict.get("generate_template", False):
        logger.warning ("Generate template files in working directory")
        generate_template (
            workflow_dir=WORKFLOW_DIR,
            templates=args_dict["generate_template"],
            workflow=args_dict["subcommand"],
            workdir=args_dict["workdir"],
            overwrite=args_dict["overwrite_template"],
            verbose=args_dict["verbose"],
            quiet=args_dict["quiet"])
        sys.exit()

    # Cluster stuff to simplify options
    if args_dict["cluster_config"]:
        logger.warning ("INITIALISING WORKFLOW IN CLUSTER MODE")
        cluster_config_fn = args_dict["cluster_config"]
        args_dict["local_cores"] = get_yaml_val(yaml_fn=cluster_config_fn, val_name="cluster_cores", default=args_dict["cores"])
        args_dict["nodes"] = get_yaml_val(yaml_fn=cluster_config_fn, val_name="cluster_nodes", default=args_dict["nodes"])
        args_dict["cluster"] = get_yaml_val(yaml_fn=cluster_config_fn, val_name="cluster_cmd", default=args_dict["cluster"])
        args_dict["config"] = args_dict["cluster_config"]
        logger.debug ("Cores:{} / Nodes:{} / Cluster_cmd:{}".format(args_dict['local_cores'], args_dict['nodes'], args_dict['cluster']))
    elif args_dict["config"] :
        logger.warning ("INITIALISING WORKFLOW IN LOCAL MODE")
    else:
        logger.error("A configuration file `--config` or a cluster configuration file `--cluster_config` is required")
        sys.exit()

    # Get and check config files
    logger.warning ("LOADING CONFIGURATIONS INFO")
    snakefile = get_snakefile_fn(workflow_dir=WORKFLOW_DIR, workflow=args_dict["subcommand"])
    configfile = get_config_fn(config=args_dict["config"])
    kwargs = filter_out_options (args_dict)
    logger.debug (kwargs)

    # Run Snakemake API
    try:
        snakemake (
            snakefile=snakefile,
            configfiles=[configfile],
            use_conda=True,
            wrapper_prefix=WRAPPER_PREFIX,
            **kwargs)
    except TypeError as E:
        logger.error("Unsupported Option Error. {}".format(E))
        sys.exit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TEST SUBPARSER FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def test_wrappers (args_dict):
    """"""
    # Cleanup data and leave
    if args_dict["clean_output"]:
        logger.info("Removing output data")
        for wrapper_name in args_dict["wrappers"]:
            wrapper_workdir = os.path.join(args_dict["workdir"], wrapper_name)
            shutil.rmtree(wrapper_workdir, ignore_errors=True)
        sys.exit()

    # Test wrappers
    for wrapper_name in args_dict["wrappers"]:
        logger.warning("Testing Wrapper {}".format(wrapper_name))
        try:
            snakefile = get_snakefile_fn(workflow_dir=WRAPPER_DIR, workflow=wrapper_name)
            wrapper_workdir = os.path.join(args_dict["workdir"], wrapper_name)
            logger.debug("Working in directory: {}".format(wrapper_workdir))

            #Run Snakemake through the API
            snakemake (
                snakefile = snakefile,
                workdir = wrapper_workdir,
                config = {"data_dir":DATA_DIR},
                wrapper_prefix = WRAPPER_PREFIX,
                use_conda = True,
                cores = args_dict["cores"],
                verbose = args_dict["verbose"],
                quiet = args_dict["quiet"])

        finally:
            logger.debug("List of file generated: {}".format(os.listdir(wrapper_workdir)))
            shutil.rmtree(os.path.join(wrapper_workdir, ".snakemake"), ignore_errors=True)
            if not args_dict["keep_output"]:
                logger.debug("Removing temporary directory")
                shutil.rmtree(wrapper_workdir, ignore_errors=True)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SCRIPT ENTRY POINT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == "__main__":
    main ()
