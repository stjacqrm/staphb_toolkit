#!/usr/bin/env python3

#authors:
# Kelsey Florek (kelsey.florek@slh.wisc.edu)
# Kevin Libuit (kevin.libuit@dgs.virginia.gov)

import sys,os,re
import argparse
from shutil import which
import shlex
import subprocess

#Class to allow printing over the output generated by nextflow
#usage:
#reprinter = Reprinter()
#reprinter.reprint("block of text")
class Reprinter:
    def __init__(self):
        self.text = ''

    def moveup(self, lines):
        for _ in range(lines):
            sys.stdout.write("\x1b[A")

    def reprint(self, text):
        # Clear previous text by overwritig non-spaces with spaces
        self.moveup(self.text.count("\n"))
        sys.stdout.write(re.sub(r"[^\s]", " ", self.text))

        # Print new text
        lines = min(self.text.count("\n"), text.count("\n"))
        self.moveup(lines)
        sys.stdout.write(text)
        self.text = text

def nextflow_printer(running_process):
    reprinter = Reprinter()
    block = ""
    running = False
    post_text = ""
    for line in iter(running_process.stdout.readline, b''):
        s = line.decode("utf-8","ignore").strip('\n')
        if "Launching" in s:
            running = True
            print(s)
        elif "process" in s:
            block += s + "\n"
        elif "executor" in s:
            continue
        elif not s:
            reprinter.reprint(block)
            block = ""
        elif s and running == False:
            print(s)
        else:
            if s == '\n':
                continue
            else:
                post_text += s + '\n'
    print(post_text)

def main():
    #get nextflow executable
    lib_path = os.path.abspath(os.path.dirname(__file__) + '/' + 'lib')
    workflows_path = os.path.abspath(os.path.dirname(__file__) + '/' + 'workflows')
    nextflow_path = os.path.join(lib_path,'nextflow')

    #setup argparser to display help if no arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(usage="staphb-wf [optional arguments] <workflow> [workflow arguments]")
    subparsers = parser.add_subparsers(title='workflows',metavar='',dest="subparser_name")

    #check if we are using docker or singularity
    if which('docker'):
        profile = 'docker'
    elif which('singularity'):
        profile = 'singularity'
    else:
        print('Singularity or Docker is not installed or not in found in PATH.')
        sys.exit(1)

    #parser for workflows
    #tredegar-----------------------------------------
    parser_tredegar = subparsers.add_parser('tredegar', help='Quality control of WGS read data.', add_help=False)
    parser_tredegar.add_argument('reads_path', type=str,help="path to the location of the reads in a fastq format")
    parser_tredegar.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"tredegar_results\".",default="tredegar_results")
    parser_tredegar.add_argument('--profile',metavar='profile_name', type=str,help="Custom nextflow profile.")

    #monroe-----------------------------------------
    parser_monroe = subparsers.add_parser('monroe', help='Consensus assembly for SARS-CoV-2 from ARTIC + Illumina protocols.', add_help=False)
    parser_monroe.add_argument('reads_path', type=str,help="path to the location of the reads in a fastq format")
    parser_monroe.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"monroe_results\".",default="monroe_results")
    parser_monroe.add_argument('--primers', type=str,choices=["V1", "V2", "V3"], help="indicate which ARTIC primers were used (V1, V2, or V3)")
    parser_monroe.add_argument('--profile',metavar='profile_name', type=str,help="Custom nextflow profile.")

    #foushee-----------------------------------------
    # parser_foushee = subparsers.add_parser('foushee', help='Reference-free SNP analysis of GAS isolates.', add_help=False)
    # parser_foushee.add_argument('reads_path', type=str,help="path to the location of the reads in a fastq format")
    # parser_foushee.add_argument('--output','-o',metavar='output', type=str,help="output directory - defaults to working directory")
    # parser_foushee.add_argument('-c',metavar='config', type=str,help="path to configuration file")
    # parser_foushee.add_argument('-t',metavar='cpus', type=int,help="number of cpus to use, defaults to 8",default=8)
    # parser_foushee.add_argument('-m',metavar='memory GB', type=int,help="number of GB of memory to use, defaults to 16",default=16)

    #dryad-----------------------------------------
    parser_dryad = subparsers.add_parser('dryad', help='A comprehensive tree building program.', add_help=False)
    parser_dryad.add_argument('reads_path', type=str,help="Path to the location of the raw reads in the fastq format.")
    parser_dryad.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"dryad_results\".",default="dryad_results")
    parser_dryad.add_argument('--core-genome','-cg',default=False, action="store_true", help="Construct a core-genome tree.")
    parser_dryad.add_argument('--snp','-s',default=False, action="store_true", help="Construct a SNP tree. Note: Requires a reference genome in fasta format (-r).")
    parser_dryad.add_argument('-ar',default=False, action="store_true", help="Detect AR mechanisms.")
    parser_dryad.add_argument('-r',metavar='<path>', type=str,help="Reference genome for SNP pipeline.")
    parser_dryad.add_argument('--profile',metavar='profile_name', type=str,help="Custom nextflow profile.")
    parser_dryad.add_argument('--sep',metavar="sep_chars",type=str,help="Dryad identifies sample names from the name of the read file by splitting the name on the specified separating characters, default \"_\".",default="_")

    args = parser.parse_args()
    if args.subparser_name == None:
        parser.print_help()

    program = args.subparser_name

    #Program specific execution code
    #-----------------------------------------
    if program == 'tredegar':
        #tredegar path
        tredegar_path = os.path.join(workflows_path,"tredegar/tredegar.nf")

        #check for user profile
        if args.profile:
            profile = args.profile

        #set work dir into local logs dir if profile not aws
        work = ""
        if profile != "aws":
            work = f"-w {args.output}/logs/work"

        #build command
        command = nextflow_path
        command = command + f" {tredegar_path} -profile {profile} -resume --reads {args.reads_path} --outdir {args.output} -with-trace {args.output}/logs/Tredegar_trace.txt -with-report {args.output}/logs/Tredegar_execution_report.html {work}"

        #run command using nextflow in a subprocess
        print("Starting the Tredegar pipeline:")
        try:
            process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
            nextflow_printer(process)
        except KeyboardInterrupt:
            print("Quitting Tredegar...")
            process.terminate()
            print("Done.")
            sys.exit(1)

    if program == 'monroe':
        #monroe path
        monroe_path = os.path.join(workflows_path,"monroe/monroe.nf")

        #check for user profile
        if args.profile:
            profile = args.profile

        #set work dir into local logs dir if profile not aws
        work = ""
        if profile != "aws":
            work = f"-w {args.output}/logs/work"

        #build command
        command = nextflow_path
        command = command + f" {monroe_path} -profile {profile} -resume --reads {args.reads_path} --primers {args.primers} --outdir {args.output} -with-trace {args.output}/logs/Monroe_trace.txt -with-report {args.output}/logs/Monroe_execution_report.html {work}"
        #run command using nextflow in a subprocess
        print("Starting the Monroe pipeline:")
        try:
            process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
            nextflow_printer(process)
        except KeyboardInterrupt:
            print("Quitting Monroe...")
            process.terminate()
            print("Done.")
            sys.exit(1)


    if program == 'foushee':
        pass

    if program == 'dryad':
        #dryad path
        dryad_path = os.path.join(workflows_path,"dryad/dryad.nf")

        #check for reference sequence
        if args.snp and args.r == None:
            parser_dryad.print_help()
            print("Please specify a reference sequence for the SNP pipeline.")
            sys.exit(1)

        #check for user profile
        if args.profile:
            profile = args.profile

        #build nextflow command
        selections = ""
        if args.ar:
            selections += " --ar"
        if args.core_genome:
            selections += " --cg"
        if args.snp:
            selections += f" --snp --snp_reference {args.r}"
        #add other arguments
        other_args = f"--name_split_on {args.sep} --outdir {args.output}"
        #build command
        command = nextflow_path
        command = command + f" {dryad_path} -profile {profile} -resume --reads {args.reads_path} {selections} {other_args}"

        #run command using nextflow in a subprocess
        print("Starting the Dryad pipeline:")
        try:
            process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
            nextflow_printer(process)
        except KeyboardInterrupt:
            print("Quitting dryad...")
            process.terminate()
            print("Done.")
            sys.exit(1)
