# -*- coding: utf-8 -*-

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Imports~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Std lib
from os.path import join

# Third party lib
import pandas as pd

# Local imports
from NanoSnake.common import *
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~check config file version~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Minimum snakemake version
config_version=1
if not "config_version" in config or config["config_version"]!= config_version:
    raise NanoSnakeError ("Wrong configuration file version. Please regenerate config with `--generate_template config -o`")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define samples sheet reference and getters~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
sample_df=pd.read_csv (config["sample_sheet"], comment="#", skip_blank_lines=True, sep="\t", index_col=0)
sample_list=list(sample_df.index)

def get_methylation_calls (wildcards):
    return sample_df.loc[wildcards.sample, "methylation_calls"]

# shortcuts to interval mode and bed file
interval_mode = config["interval_mode"]
external_bed = config["external_bed"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define IO for each rule~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Build output files dictionnary
input_d=defaultdict(OrderedDict)
output_d=defaultdict(OrderedDict)
log_d=OrderedDict()

rule_name="get_genome"
if config["genome"].startswith("ftp"):
    input_d[rule_name]["ref"]=FTP.remote(config["genome"])
elif config["genome"].startswith("http"):
    input_d[rule_name]["ref"]=HTTP.remote(config["genome"])
else:
    input_d[rule_name]["ref"]=config["genome"]
output_d[rule_name]["ref"]=join("results","input","genome","ref.fa")
output_d[rule_name]["index"]=join("results","input","genome","ref.fa.fai")
log_d[rule_name]=join("logs",rule_name,"ref.log")

if interval_mode == "cpg_islands":
    rule_name="pycometh_cgi_finder"
    input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
    output_d[rule_name]["tsv"]=join("results","methylation","pycometh_cgi_finder","CGI.tsv.gz")
    output_d[rule_name]["bed"]=join("results","methylation","pycometh_cgi_finder","CGI.bed")
    output_d[rule_name]["bed_index"]=join("results","methylation","pycometh_cgi_finder","CGI.bed.idx")
    log_d[rule_name]=join("logs",rule_name,"ref.log")

rule_name="pycometh_cpg_aggregate"
input_d[rule_name]["tsv"]=get_methylation_calls
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
output_d[rule_name]["tsv"]=join("results","methylation","pycometh_cpg_aggregate","{sample}.tsv.gz")
output_d[rule_name]["bed"]=join("results","methylation","pycometh_cpg_aggregate","{sample}.bed")
output_d[rule_name]["bed_index"]=join("results","methylation","pycometh_cpg_aggregate","{sample}.bed.idx")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="pycometh_interval_aggregate"
input_d[rule_name]["tsv"]=output_d["pycometh_cpg_aggregate"]["tsv"]
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
# Conditional input annotation file
if interval_mode == "cpg_islands":
    input_d[rule_name]["annot"]=output_d["pycometh_cgi_finder"]["bed"]
elif interval_mode == "cpg_islands" and external_bed:
    input_d[rule_name]["annot"]=external_bed
output_d[rule_name]["tsv"]=join("results","methylation","pycometh_interval_aggregate","{sample}.tsv.gz")
output_d[rule_name]["bed"]=join("results","methylation","pycometh_interval_aggregate","{sample}.bed")
output_d[rule_name]["bed_index"]=join("results","methylation","pycometh_interval_aggregate","{sample}.bed.idx")
log_d[rule_name]=join("logs",rule_name,"{sample}.log")

rule_name="pycometh_meth_comp"
input_d[rule_name]["tsv"]=expand(join("results","methylation","pycometh_interval_aggregate","{sample}.tsv.gz"), sample=sample_list) # Aggregate samples
input_d[rule_name]["ref"]=output_d["get_genome"]["ref"]
output_d[rule_name]["tsv"]=join("results","methylation","pycometh_meth_comp","meth_comp.tsv.gz")
output_d[rule_name]["bed"]=join("results","methylation","pycometh_meth_comp","meth_comp.bed")
output_d[rule_name]["bed_index"]=join("results","methylation","pycometh_meth_comp","meth_comp.bed.idx")
log_d[rule_name]=join("logs",rule_name,"meth_comp.log")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define all output depending on config file~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

all_output=[]
for rule_name, rule_d in output_d.items():
    for option, output in rule_d.items():
        all_output.append(output)
all_output = flatten_list(all_output)

all_expand = []
for output in all_output:
    all_expand.append(list(set(expand(output, sample=sample_list))))
all_expand = flatten_list(all_expand)

print (all_expand)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RULES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rule all:
    input: all_expand

rule_name="get_genome"
rule get_genome:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "get_genome"

rule_name="pycometh_cgi_finder"
rule pycometh_cgi_finder:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_cgi_finder"

rule_name="pycometh_cpg_aggregate"
rule pycometh_cpg_aggregate:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_cpg_aggregate"

rule_name="pycometh_interval_aggregate"
rule pycometh_interval_aggregate:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_interval_aggregate"

rule_name="pycometh_meth_comp"
rule pycometh_meth_comp:
    input: **input_d[rule_name]
    output: **output_d[rule_name]
    log: log_d[rule_name]
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    wrapper: "pycometh_meth_comp"
