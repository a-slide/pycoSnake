# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
import sys
from collections import *
import shutil
import yaml
import inspect

# Third party lib
import pandas as pd
from loguru import logger as log

#~~~~~~~~~~~~~~CUSTOM EXCEPTION CLASS~~~~~~~~~~~~~~#
class NanoSnakeError (Exception):
    """ Basic exception class"""
    pass

class NanoSnakeWarning (Warning):
    """ Basic Warning class"""
    pass

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#

def mkdir (fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    try:
        os.makedirs (fn, exist_ok=exist_ok)
    except:
        raise NanocomporeError ("Error creating output folder `{}`".format(fn))

def access_file (fn, **kwargs):
    """ Check if the file is readable """
    return os.path.isfile (fn) and os.access (fn, os.R_OK)

def get_log_level(quiet=False, verbose=False):
    """return logging level depending on verbosity args"""
    if quiet:
        return "WARNING"
    if verbose:
        return "DEBUG"
    return "INFO"

def set_log_level(quiet=False, verbose=False):
    level = get_log_level(quiet,verbose)
    log.remove()
    log.add (sys.stderr, level=level)

def jhelp (f:"python function or method"):
    """
    Display a Markdown pretty help message for functions and class methods (default __init__ is a class is passed)
    jhelp also display default values and type annotations if available.
    Undocumented options are not displayed.
    The docstring synthax should follow the markdown formated convention below
    * f
        Function or method to display the help message for
    """
    # For some reason signature is not always importable. In these cases the build-in help is called instead
    try:
        from IPython.core.display import display, Markdown, HTML
        import inspect
    except (NameError, ImportError) as E:
        NanocomporeWarning ("jupyter notebook is required to use this function. Please verify your dependencies")
        help(f)
        return

    if inspect.isclass(f):
        f = f.__init__

    if inspect.isfunction(f) or inspect.ismethod(f):

        # Parse arguments default values and annotations
        sig_dict = OrderedDict()
        for name, p in inspect.signature(f).parameters.items():
            sig_dict[p.name] = []
            # Get Annotation
            if p.annotation != inspect._empty:
                sig_dict[p.name].append(": {}".format(p.annotation))
            # Get default value if available
            if p.default == inspect._empty:
                sig_dict[p.name].append("(required)")
            else:
                sig_dict[p.name].append("(default = {})".format(p.default))

        # Parse the docstring
        doc_dict = OrderedDict()
        descr = []
        lab=None
        for l in inspect.getdoc(f).split("\n"):
            l = l.strip()
            if l:
                if l.startswith("*"):
                    lab = l[1:].strip()
                    doc_dict[lab] = []
                elif lab:
                    doc_dict[lab].append(l)
                else:
                    descr.append(l)

        # Reformat collected information in Markdown synthax
        s = "---\n\n**{}.{}**\n\n{}\n\n---\n\n".format(f.__module__, f.__name__, " ".join(descr))
        for k, v in doc_dict.items():
            s+="* **{}** *{}*\n\n{}\n\n".format(k, " ".join(sig_dict[k]), " ".join(v))

        # Display in Jupyter
        display (Markdown(s))

#~~~~~~~~~~~~~~SNAKEMAKE HELPER FUNCTIONS~~~~~~~~~~~~~~#

def all_in (collection, val_list):
    """"""
    for val in val_list:
        if not val in collection:
            return False
    return True

def get_threads (config, rule_name, default=1):
    try:
        return config[rule_name]["threads"]
    except (KeyError, TypeError):
        return default

def get_opt (config, rule_name, default=""):
    try:
        return config[rule_name]["opt"]
    except KeyError:
        return default

def get_mem (config, rule_name, default=1000):
    try:
        return config[rule_name]["mem"]
    except KeyError:
        return default

def get_output (out, rule_name, val):
    try:
        return out[rule_name][val]
    except KeyError:
        return val

#~~~~~~~~~~~~~~MAIN HELPER FUNCTIONS~~~~~~~~~~~~~~#

def add_argument_group (parser, title):
    """Add group only is it doesn't exist yet"""
    for group in parser._action_groups:
        if group.title == title:
            return group
    return parser.add_argument_group(title)

def unlock_dir (workdir):
    """"""
    lockdir = os.path.join(workdir, ".snakemake", "locks")
    try:
        shutil.rmtree(lockdir)
    except:
        raise NanoSnakeError ("Cannot remove lock")

def generate_template (workflow_dir, templates, workflow, workdir="./", overwrite=False, verbose=False, quiet=False):
    """"""
    set_log_level(quiet=quiet, verbose=verbose)

    templates_to_fname = {
        "sample_sheet":"sample_sheet.tsv" ,
        "config":"config.yaml",
        "cluster_config":"cluster_config.yaml"}

    for template, fname in templates_to_fname.items():
        if template in templates or "all" in templates:

            # Create src path and test readability
            src_fn = os.path.join (workflow_dir, workflow, "templates", fname)

            # Create destination file name and test if it exists
            dest_fn = os.path.join(workdir, fname)
            if not os.path.isfile(dest_fn) or overwrite:
                shutil.copy2 (src_fn, dest_fn)
                log.info(f"{template} file created in working directory ")
            else:
                log.info(f"{template} file already exist in working directory. Use --overwrite_template to overwrite existing file")

def get_yaml_val(yaml_fn, val_name, default):
    """"""
    try:
        with open(yaml_fn) as fp:
            y = yaml.load(fp, Loader=yaml.FullLoader)
            return y[val_name]
    except:
        return default

def get_config_fn (config):
    """"""
    # Try loading the config file
    try:
        with open(config) as fp:
            yaml.load(fp, Loader=yaml.FullLoader)
    except:
        raise NanoSnakeError ("The provided config file is not readeable or not a valid yaml file")

    return os.path.abspath(config)

def get_snakefile_fn (workflow_dir, workflow):
    """"""
    snakefile = os.path.join (workflow_dir, workflow, "snakefile.py")
    if not access_file(snakefile):
        raise NanoSnakeError ("The snakefile file is not readable")
    return os.path.abspath(snakefile)

def get_sample_sheet (sample_sheet, required_fields=[]):
    """"""
    if not sample_sheet:
        raise NanoSnakeError ("A sample_sheet file (--sample_sheet) is required to run the workflow")
    try:
        sample_df = pd.read_csv (sample_sheet, comment="#", skip_blank_lines=True, sep="\t")
    except:
        raise NanoSnakeError ("Cannot open provided sample sheet")
    for f in required_fields:
        if not f in sample_df.columns:
            raise NanoSnakeError ("The provided sample sheet does not contain the required fieds: {}".format(" ".join(required_fields)))
    return os.path.abspath(sample_sheet)

def required_option (name, var):
    """"""
    if not var:
        raise NanoSnakeError (f"Option --{name} is required to run the workflow")
    return var
