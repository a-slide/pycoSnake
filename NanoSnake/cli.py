#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
import os
import sys
import argparse
import pkg_resources
import logging
import yaml
import shutil
from collections import *

# Third party library
from snakemake import snakemake

# Local imports
from NanoSnake import __version__ as package_version
from NanoSnake import __name__ as package_name
from NanoSnake import __description__ as package_description

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~LOGGING INFO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)
logLevel_dict = {2:logging.DEBUG, 1:logging.INFO, 0:logging.WARNING}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GLOBAL DIRS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ENVS_DIR = pkg_resources.resource_filename (package_name, 'envs')
SCRIPTS_DIR = pkg_resources.resource_filename (package_name, 'scripts')
RULES_DIR = pkg_resources.resource_filename (package_name, 'rules')
SNAKEFILES_DIR = pkg_resources.resource_filename (package_name, 'snakefiles')
TEMPLATES_DIR = pkg_resources.resource_filename (package_name, 'templates')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CLI ENTRY POINT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def main(args=None):
    """
    Main entry point for NanoSnake command line interface
    """
    # Parse common options first
    common_parser = argparse.ArgumentParser (add_help=False)
    common_parser_IO = common_parser.add_argument_group('Common input/output options')
    common_parser_IO.add_argument("--sample_sheet", "-s", required=True, type=str,
        help="Path to a tabulated sample sheet (--sample_template to generate a template).\
        Required if not given in the configuration file")
    common_parser_IO.add_argument("--config_file", "-c", default="", type=str,
        help="Snakemake configuration YAML file (--config_template to generate a template).\
        If not given the default package configuration file is used")
    common_parser_IO.add_argument("--working_dir", "-d", default="./", type=str,
        help="Path to the working dir where to deploy the worflow (default: %(default)s)")
    common_parser_template = common_parser.add_argument_group('Generate template files')
    common_parser_template.add_argument("--config_template", action="store_true",
        help="Output a template YAML configuration file in the working directory and quit)")
    common_parser_template.add_argument("--sample_template", action="store_true",
        help="Output a template TSV sample sheet file in the working directory and quit)")
    common_parser_other = common_parser.add_argument_group('Other options')
    common_parser_other.add_argument("--verbose_level", "-v", type=int, default=1, choices=[0,1,2],
        help="Level of verbosity, from 2 (Chatty) to 0 (Nothing) (default: %(default)s)")

    # Parser and subparsers for command
    parser = argparse.ArgumentParser ()
    parser.add_argument('--version', action='version', version="{} v{}".format(package_name, package_version))
    subparsers = parser.add_subparsers (description='NanoSnake implements the following subcommands', dest='subcommands')
    subparsers.required = True

    # DNA_methylation subparser
    subparser_dm = subparsers.add_parser('DNA_methylation', parents=[common_parser],
        description="Workflow to evaluate DNA methylation in Nanopore data using Nanopolish")
    # Additional optional agrs that can also be provided in the config file
    subparser_dm.add_argument("--reference", "-r", required=True, type=str,
        help="Path to a Fasta reference file to be used for read mapping.\
        Required if not given in the configuration file)")

    # RNA_count subparser
    subparser_rc = subparsers.add_parser('RNA_count', parents=[common_parser],
        description="Workflow aligning RNA on the transcriptome an generating estimated counts")

    # Parse args and call subcommands
    args = parser.parse_args()
    args_dict = {i:j for i,j in  vars(args).items() if i != "subcommands"}
    if args.subcommands == "DNA_methylation":
        DNA_methylation (**args_dict)
    elif args.subcommands == "RNA_count":
        RNA_count (**args_dict)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SUBPARSERS FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def DNA_methylation (sample_sheet, reference, config_file="", working_dir="./", config_template="", sample_template="", verbose_level=1):

    # Define loglevel
    logger.setLevel (logLevel_dict.get (verbose_level, logging.INFO))

    # Define default package files and dir
    snake_fn = os.path.join (SNAKEFILES_DIR, 'DNA_methylation_snakefile')
    default_config = os.path.join (TEMPLATES_DIR, 'DNA_methylation_config.yaml')
    default_sample_sheet = os.path.join (TEMPLATES_DIR, 'DNA_methylation_samples.tsv')

    # If templates required
    if config_template:
        file_dest = os.path.join(working_dir, "config.yaml")
        logger.info (f"Create template config file: {file_dest}")
        shutil.copyfile (default_config, file_dest)
    elif sample_template:
        file_dest = os.path.join(working_dir, "samples.tsv")
        logger.info (f"Create template config file: {file_dest}")
        shutil.copyfile (default_sample_sheet, file_dest)

    # Load the config file and run snakemake through the API
    else:
        # Add command line args to config or get if from the config yaml file
        logger.info ("Build config dict for snakemake")
        config = {
            "sample_sheet": sample_sheet,
            "reference": reference,
            "envs_dir": ENVS_DIR,
            "scripts_dir": SCRIPTS_DIR
            "rules_dir": RULES_DIR}
        logger.debug (config)

        # Run Snakemake through the API
        logger.warning ("RUN SNAKEMAKE PIPELINE")
        config_file =
        snakemake (snakefile=snake_fn, configfile=config_file, config=config, workdir=working_dir)


def RNA_count (args):
    pass

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SCRIPT ENTRY POINT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if __name__ == "__main__":
    main ()
