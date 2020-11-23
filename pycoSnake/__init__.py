# -*- coding: utf-8 -*-

__name__ = "pycoSnake"
__version__ = "0.2.4.dev3"
__description__ = "pycoSnake is a neatly wrapped collection of snakemake workflows for analysing Illumina and nanopore sequencing datasets. It is easy to install with conda and simple to run on a local computer or in a cluster environment"
__url__ = "https://github.com/a-slide/pycoSnake"
__licence__ = "MIT"
__author__ = "Adrien Leger"

workflows_info = {
    "DNA_ONT" : {
        "version" : "0.4",
        "description" : "Workflow for DNA analysis of Nanopore data including, alignment, methylation calling and SV calling"},
    "RNA_illumina" : {
        "version" : "0.1",
        "description" : "Workflow for RNA analysis of Illunina data including, alignment and gene abundance estimation"},
    "test_wrappers" : {
        "version" : "0.4",
        "description" : "Test Nanosnake wrappers"}}
