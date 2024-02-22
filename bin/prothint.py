#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2019, Georgia Institute of Technology, USA
#
# ProtHint: Pipeline for generating genome wide footprints of homologous
# proteins
# ==============================================================

import argparse
import os
import sys
import subprocess
import time
import shutil
import tempfile


VERSION = '2.6.0'


def stderr_message(msg):
    sys.stderr.write(f"[{time.ctime()}] {msg}\n")


def main():
    args = parseCmd()

    args = setEnvironment(args)

    if not args.diamond_db:
        args = processInputProteins(args)

    if args.geneSeeds and args.prevGeneSeeds:
        nextIteration(args)
    else:
        standardRun(args)


def standardRun(args):
    """Execute a standard ProtHint run

    Args:
        args: Command line arguments
    """
    if not args.geneSeeds:
        args.geneSeeds = runGeneMarkES(args)
    else:
        stderr_message('''Skipping GeneMark-ES, using the supplied gene seeds
            file instead''')

    translateSeeds(args, args.geneSeeds)

    if not args.diamondPairs:
        args.diamondPairs = runDiamond(args)
    else:
        stderr_message('''Skipping DIAMOND, using the supplied DIAMOND output
            file instead''')

    prepareSeedSequences(args)

    runSpaln(args, args.diamondPairs, args.pbs, args.minExonScore,
             args.nonCanonicalSpliceSites, args.longGene, args.longProtein)

    checkOutputs(args)
    flagTopProteins(args)
    processSpalnOutput(args)

    if args.cleanup:
        cleanup(args)
    if not args.diamond_db:
        os.remove(args.proteins)  # Args contains not the original file, but temporary
    stderr_message('ProtHint finished.')


def nextIteration(args):
    """Run a next iteration of ProtHint. ProtHint is only run for gene
    seeds which are new or modified in the --geneSeeds file compared to
    --prevGeneSeeds. Hints for genes which are the same are reused from the
    --prevSpalnGff file but their gene_ids are updated to match the new
    seed file.

    Args:
        args: Command line arguments
    """
    stderr_message("ProtHint is running in the iterative mode.")
    prepareDataForNextIteration(args)

    diamondPairs = ""
    if os.path.getsize("uniqueSeeds.gtf") != 0:
        translateSeeds(args, "uniqueSeeds.gtf")
        diamondPairs = runDiamond(args)
        prepareSeedSequences(diamondPairs)
        runSpaln(diamondPairs, args.pbs, args.minExonScore,
                 args.nonCanonicalSpliceSites, args.longGene, args.longProtein)
        flagTopProteins(diamondPairs)
        # Append subset of hints from the previous iteration to the current result
        os.chdir(args.workdir)
        with open("Spaln/spaln.gff", "a") as new:
            with open("prevHints.gff", "r") as prev:
                for line in prev:
                    new.write(line)
    else:
        sys.stderr.write("Warning: No unique gene seeds were detected in the " +
                         args.geneSeeds + " input file. ProtHint will only " \
                         "update seed gene IDs of hints to match the IDs in " \
                         "the new seed gene file.\n")
        if not os.path.isdir("Spaln"):
            os.mkdir("Spaln")
        shutil.move("prevHints.gff", "Spaln/spaln.gff")

    processSpalnOutput(args)

    os.remove(args.proteins)  # Args contains not the original file, but temporary
    sys.stderr.write("[" + time.ctime() + "] ProtHint finished.\n")


def prepareDataForNextIteration(args):
    """Select gene seeds which are unique in this iteration and hints from
    previous iteration of ProtHint which correspond to seeds which are
    identical in both files (with gene ids updated to match the ids in the
    new seed file.

    Args:
        args.geneSeeds (filepath): Path to current gene seeds
        args.prevGeneSeeds (filepath): Path to previous genes seeds
        args.prevSpalnGff (filepath): Path to previous scored hints
    """
    stderr_message("Selecting a subset of data to run in the next iteration")

    os.chdir(args.workdir)

    callScript('print_longest_isoform.py -i', f'''{args.geneSeeds}
        -o longest_seed_isoforms.gtf''')

    callScript('print_longest_isoform.py -i', f'''{args.prevGeneSeeds}
        -o longest_prevSeed_isoforms.gtf''')

    callScript("select_for_next_iteration.py",
        ' '.join(["--geneSeeds ", "longest_seed_isoforms.gtf",
            "--prevGeneSeeds", "longest_prevSeed_isoforms.gtf",
            "--prevSpalnGff", args.prevSpalnGff,
            "--uniqueNewSeedsOut", "uniqueSeeds.gtf"
            "--identicalSpalnGffOut", "prevHints.gff"])
        )
    os.remove("longest_prevSeed_isoforms.gtf")
    os.remove("longest_seed_isoforms.gtf")


def runGeneMarkES(args):
    """Run GeneMark-ES

    Args:
        pbs (boolean): Whether to run on pbs
        fungus (boolean): Whether to run in fungus mode

    Returns:
        string: Path to genemark gff output
    """
    stderr_message("Running GeneMark-ES.")
    ESDir = args.workdir + "/GeneMark_ES"
    if not os.path.isdir(ESDir):
        os.mkdir(ESDir)
    os.chdir(ESDir)

    cmd = ["gmes_petap.pl", "--verbose", "--cores", args.threads,
        "--ES", "--seq", args.genome, "--soft", "auto"]

    if args.pbs:
        cmd.append('--pbs')
    if args.fungus:
        cmd.append('--fungus')

    runSubprocess(cmd)  # gmes_petap.pl should be in your PATH.
    '''callDependency("gmes_petap.pl", "--verbose --cores " + threads + pbsFlag +
                   " --ES --seq " + genome + " --soft auto" + fungusFlag,
                   "GeneMarkES")'''

    stderr_message("GeneMark-ES finished.")

    dir_list = ['data', 'info', 'run', 'output']
    for folder in dir_list:
        if os.path.isdir(folder):
            shutil.rmtree(folder)

    return os.path.abspath("genemark.gtf")


def translateSeeds(args, seedGenes):
    """Translate GeneMark-ES seeds to proteins

    Args:
        seedGenes (filepath): Path to file with seed genes
    """
    stderr_message("Translating gene seeds to proteins")
    os.chdir(args.workdir)

    callScript('print_longest_isoform.py', f'-i {seedGenes}')

    systemCall("grep \tCDS\t longest_seed_isoforms.gtf > " +
               "longest_seed_isoforms_cds.gtf")
    os.remove("longest_seed_isoforms.gtf")

    callScript("proteins_from_gtf.pl", "--stat gene_stat.yaml --seq " +
               args.genome + " --annot longest_seed_isoforms_cds.gtf --out " +
               "seed_proteins.faa --format GTF")
    os.remove("longest_seed_isoforms_cds.gtf")

    stderr_message("Translation of seeds finished")


def runDiamond(args):
    """Run DIAMOND protein search

    Args:
        maxProteins (int): Maximum number of protein hits per seed gene.
        evalue (float): Maximum e-value of DIAMOND hits

    Returns:
        string: Path to DIAMOND output
    """
    stderr_message("Running DIAMOND")

    diamondDir = args.workdir + "/diamond"
    if not os.path.isdir(diamondDir):
        os.mkdir(diamondDir)
    os.chdir(diamondDir)

    if not args.diamond_db:
        # Make DIAMOND db
        callDependency("diamond", "makedb --in " + args.proteins + " -d diamond_db " +
                       "--threads " + args.threads)
        args.diamond_db = 'diamond_db'

    # Actual DIAMOND run
    callDependency("diamond", "blastp --query ../seed_proteins.faa --db " +
                   f"{args.diamond_db} --outfmt 6 qseqid sseqid --out diamond.out " +
                   "--max-target-seqs " + str(args.maxProteinsPerSeed) + " --max-hsps 1 " +
                   "--threads " + args.threads + " --evalue " + str(args.evalue))

    stderr_message('DIAMOND finished')
    return os.path.abspath("diamond.out")


def prepareSeedSequences(args):
    """Prepare nucleotide sequences for seed genes

    Args:
        diamondPairs (filepath): Path to file with seed gene-protein pairs
    """
    sys.stderr.write("[" + time.ctime() + "] Preparing pairs for alignments\n")
    os.chdir(args.workdir)

    callScript("nucseq_for_selected_genes.pl", "--seq " + args.genome + " --out " +
               "nuc.fasta --gene gene_stat.yaml --list " + args.diamondPairs)

    stderr_message('Preparation of pairs finished')


def runSpaln(args, diamondPairs, pbs, minExonScore, nonCanonical,
             longGene, longProtein):
    """Run Spaln spliced alignment and score the outputs with spaln-boundary-scorer

    Args:
        diamondPairs (filePath): Path to file with seed gene-protein pairs
                                 to align
        pbs (boolean): Whether to run on pbs
        minExonScore (float): Discard all hints inside/neighboring exons with
                              score lower than minExonScore
        nonCanonical (bool): Whether to predict non-canonical introns
        longGene (int): Threshold for what is considered a long gene in
                        Spaln alignment
        longProtein (int): Threshold for what is considered a long protein in
                           Spaln alignment
    """
    spalnDir = os.path.join(args.workdir, "Spaln")
    if not os.path.isdir(spalnDir):
        os.mkdir(spalnDir)
    os.chdir(spalnDir)

    nonCanonicalFlag = ""
    if nonCanonical:
        nonCanonicalFlag = " --nonCanonical "

    if not pbs:
        callScript("run_spliced_alignment.pl", "--cores " + args.threads +
                   " --nuc ../nuc.fasta --list " + args.diamondPairs + " --prot " +
                   args.proteins + " --v --aligner spaln --min_exon_score " +
                   str(minExonScore) + nonCanonicalFlag +
                   " --longGene " + str(longGene) +
                   " --longProtein " + str(longProtein))
    else:
        callScript("run_spliced_alignment_pbs.pl", "--N 120 --K " + args.threads +
                   " --seq ../nuc.fasta --list " + args.diamondPairs + " --db " +
                   args.proteins + " --v --aligner spaln --min_exon_score " +
                   str(minExonScore) + nonCanonicalFlag +
                   " --longGene " + str(longGene) +
                   " --longProtein " + str(longProtein))


def checkOutputs(args):
    """Check whether all intermediate outputs were correctly created

    Args:
        args.diamondPairs (filepath): Path to file with seed gene-protein pairs
    """
    os.chdir(args.workdir)

    msg = 'This error can be caused by:\n' \
          '    a) The set of input proteins is too small and/or the ' \
          'proteins are too remote.\n' \
          '    b) The gene seeds identified by GeneMark-ES (or the ' \
          'supplied gene seeds in case the option "--geneSeeds" was ' \
          'used) are incorrect. Please try running GeneMark-ES ' \
          'separately to identify errors related to gene seeds (' \
          'https://github.com/gatech-genemark/ProtHint#genemark-es).'

    if os.stat(args.diamondPairs).st_size == 0:
        sys.exit('error: No homologous proteins were found by DIAMOND (' +
                 args.diamondPairs + ' is empty).\nThis error can be caused' + msg)

    if os.stat("Spaln/spaln.gff").st_size == 0:
        sys.exit('error: No spliced alignments were created by Spaln (' +
                 args.workdir + '/Spaln/spaln.gff is empty).\n' + msg)


def flagTopProteins(args):
    """Label hints which were mapped from the best DIAMOND target

    Args:
        diamondPairs (filepath): Path to file with seed gene-protein pairs
    """
    stderr_message('Flagging top chains')
    os.chdir(args.workdir)

    callScript("flag_top_proteins.py", "Spaln/spaln.gff " + args.diamondPairs +
               " > tmp")
    shutil.move("tmp", "Spaln/spaln.gff")


def processSpalnOutput(args):
    """Prepare the final output from Spaln result scored by spaln-boundary-scorer
       Convert the output to GeneMark and Augustus compatible formats

    Args:
        nonCanonical (bool): Whether to add non-canonical introns to the
                             high-confidence set
    """
    stderr_message('Processing the output')
    os.chdir(args.workdir)

    processSpalnIntrons()
    processSpalnStops()
    processSpalnStarts()
    systemCall("sort -k1,1 -k4,4n -k5,5n -o prothint.gff prothint.gff")

    printTopChains()

    # High confidence
    nonCanonicalFlag = ""
    if args.nonCanonicalSpliceSites:
        nonCanonicalFlag = " --addAllSpliceSites "

    callScript("print_high_confidence.py", "prothint.gff " + nonCanonicalFlag +
               " > evidence.gff")

    # Augustus compatible format
    callScript("prothint2augustus.py", "prothint.gff evidence.gff "
               "top_chains.gff prothint_augustus.gff")
    stderr_message('Output processed')


def processSpalnIntrons():
    systemCall("grep Intron Spaln/spaln.gff > introns.gff || [[ $? == 1 ]]")

    # Filter out introns with alignment score < 0.1
    callScript("print_high_confidence.py", "introns.gff --intronCoverage 0 " +
               "--intronAlignment 0.1 --addAllSpliceSites > introns_01.gff")
    os.remove("introns.gff")

    callScript("combine_gff_records.pl", "--in_gff introns_01.gff --out_gff " +
               "prothint.gff")
    os.remove("introns_01.gff")


def processSpalnStops():
    systemCall("grep stop_codon Spaln/spaln.gff > stops.gff || [[ $? == 1 ]]")

    # Filter out stops with alignment score < 0.01
    callScript("print_high_confidence.py", "stops.gff --stopCoverage 0 " +
               "--stopAlignment 0.01 > stops_01.gff")
    os.remove("stops.gff")

    callScript("combine_gff_records.pl", "--in_gff stops_01.gff --out_gff " +
               "stops_01_combined.gff")
    systemCall("cat stops_01_combined.gff >> prothint.gff")

    os.remove("stops_01.gff")
    os.remove("stops_01_combined.gff")


def processSpalnStarts():
    systemCall("grep start_codon Spaln/spaln.gff > starts.gff || [[ $? == 1 ]]")

    # Filter out starts with alignment score < 0.01
    callScript("print_high_confidence.py", "starts.gff --startCoverage 0 " +
               "--startAlignment 0.01 > starts_01.gff")
    os.remove("starts.gff")

    callScript("combine_gff_records.pl", "--in_gff starts_01.gff --out_gff " +
               "starts_01_combined.gff")
    os.remove("starts_01.gff")

    # The rest of this function counts CDS overlap of starts

    systemCall("sort -k1,1 -k4,4n -k5,5n starts_01_combined.gff > " +
               "starts_01_combined_sorted.gff")
    os.remove("starts_01_combined.gff")

    systemCall("grep CDS Spaln/spaln.gff > cds.gff || [[ $? == 1 ]]")
    callScript("combine_gff_records.pl", "--in_gff cds.gff --out_gff " +
               "cds_combined.gff")
    os.remove("cds.gff")

    # Only count CDS regions which have an upstream support
    # (by start codon or intron) in hints
    callScript("cds_with_upstream_support.py", "cds_combined.gff " +
               "starts_01_combined_sorted.gff prothint.gff > tmp")
    shutil.move("tmp", "cds_combined.gff")

    systemCall("sort -k1,1 -k4,4n -k5,5n cds_combined.gff > " +
               "cds_combined_sorted.gff")
    os.remove("cds_combined.gff")

    callScript("count_cds_overlaps.py", "starts_01_combined_sorted.gff " +
               "cds_combined_sorted.gff >> prothint.gff")

    os.remove("starts_01_combined_sorted.gff")
    os.remove("cds_combined_sorted.gff")


def printTopChains():
    systemCall("grep topProt=TRUE Spaln/spaln.gff > topProteins.gff " +
               "|| [[ $? == 1 ]]")

    if os.stat("topProteins.gff").st_size == 0:
        sys.exit('error: The "topProt=TRUE" flag is missing in the '
                 'Spaln/spaln.gff output file. This issue can be caused by '
                 'the presence of special characters in the fasta headers of '
                 'input files. Please remove any special characters and '
                 're-run ProtHint. See https://github.com/gatech-genemark/ProtHint#input '
                 'for more details about the input format.')

    callScript("print_high_confidence.py", "topProteins.gff --startCoverage 0 " +
               "--startAlignment 0.01 --stopCoverage 0 --stopAlignment 0.01 " +
               "--intronCoverage 0 --intronAlignment 0.1 --addAllSpliceSites" +
               " > topProteinsFiltered.gff")
    os.remove("topProteins.gff")

    callScript("make_chains.py", "topProteinsFiltered.gff > top_chains.gff")
    os.remove("topProteinsFiltered.gff")


def cleanup(args):
    """Delete temporary files and intermediate results
    """
    sys.stderr.write("[" + time.ctime() + "] Cleaning up\n")
    os.chdir(args.workdir)

    try:
        os.remove("gene_stat.yaml")
        # os.remove("seed_proteins.faa")
        os.remove("nuc.fasta")
    except OSError:
        pass

    try:
        os.remove("diamond/diamond_db.dmnd")
    except OSError:
        pass

    try:
        os.remove("Spaln/spaln.gff")
    except OSError:
        pass


def processInputProteins(args):
    """Remove dots from the input file with proteins.
       OrhoDB protein sequences sometimes end with a dot. This format is not
       compatible with DIAMOND.
       Clean fasta headers by removing pipe ("|") characters.
    """
    stderr_message('Pre-processing protein input.')
    os.chdir(args.workdir)
    protFile = tempfile.NamedTemporaryFile(delete=False, dir='.', prefix="prot")
    systemCall('sed \"s/\.//\" ' + args.proteins + ' | sed \"s/|/_/g\" > ' +
               protFile.name)

    args.proteins = checkFileAndMakeAbsolute(protFile.name)

    return args


def setEnvironment(args):
    """Set up and check variables

    Args:
        args (dictionary): Command line arguments

    """
    ProtHintRef = 'https://doi.org/10.1093/nargab/lqaa026'
    DIAMONDRef = 'https://doi.org/10.1038/nmeth.3176'
    SpalnRef = 'https://doi.org/10.1093/bioinformatics/btn460'

    sys.stderr.write("ProtHint Version " + VERSION + "\n")
    sys.stderr.write("Copyright 2019, Georgia Institute of Technology, USA\n\n")
    sys.stderr.write("Please cite\n")
    sys.stderr.write("  - ProtHint: " + ProtHintRef + "\n")
    sys.stderr.write("  - DIAMOND:  " + DIAMONDRef + "\n")
    sys.stderr.write("  - Spaln:    " + SpalnRef + "\n\n")
    global binDir
    args.workdir = os.path.abspath(args.workdir)
    binDir = os.path.abspath(os.path.dirname(__file__))

    args.genome = checkFileAndMakeAbsolute(args.genome)
    args.proteins = checkFileAndMakeAbsolute(args.proteins)

    if args.ProSplign:
        sys.exit("error: ProSplign is not supported in this version of ProtHint. "
                 "For running ProtHint with ProSplign, use ProtHint release v2.4.0: "
                 "https://github.com/gatech-genemark/ProtHint/releases/tag/v2.4.0")

    if args.geneMarkGtf:
        if args.geneSeeds:
            sys.exit("error: please specify either --geneSeeds or\n"
                     "--geneMarkGtf. The arguments are identical,\n"
                     "--geneMarkGtf is supported for backwards compatibility.")
        args.geneSeeds = args.geneMarkGtf

    if args.geneSeeds:
        args.geneSeeds = checkFileAndMakeAbsolute(args.geneSeeds)

    if args.prevGeneSeeds:
        if not args.prevSpalnGff or not args.prevSpalnGff:
            sys.exit("error: --prevSpalnGff and --geneSeeds must be given when\n"
                     "using --prevGeneSeeds")
        args.prevGeneSeeds = checkFileAndMakeAbsolute(args.prevGeneSeeds)

    if args.prevSpalnGff:
        if not args.geneSeeds or not args.prevGeneSeeds:
            sys.exit("error: --geneSeeds and --prevGeneSeeds options must\n"
                     "be given when using --prevSpalnGff option")
        args.prevSpalnGff = checkFileAndMakeAbsolute(args.prevSpalnGff)

    if args.diamondPairs:
        if not args.geneSeeds:
            sys.exit("error: File with gene predictions (--geneSeeds\n"
                     "option) not specified. When using the --diamondPairs\n"
                     "option, a prediction file with seed genes corresponding\n"
                     "to seed genes in the DIAMOND pairs file must be specified.")
        args.diamondPairs = checkFileAndMakeAbsolute(args.diamondPairs)

    if not os.path.isdir(args.workdir):
        os.mkdir(args.workdir)

    # Log info about cmd
    callDir = "Called from: " + os.path.abspath(".") + "\n"
    cmd = "Cmd: " + " ".join(sys.argv) + "\n\n"
    sys.stderr.write(callDir)
    sys.stderr.write(cmd)
    with open(os.path.join(args.workdir, "cmd.log"), "w") as file:
        file.write(callDir)
        file.write(cmd)

    if args.threads < 1:
        if hasattr(os, "sched_getaffinity"):
            args.threads = len(os.sched_getaffinity(0))
        else:
            args.threads = os.cpu_count() or 1
    args.threads = str(args.threads)

    return args


def checkFileAndMakeAbsolute(file):
    """Check if a file with given name exists and make the path to it absolute

    Args:
        file (filepath): Path to the file
    """
    if not os.path.isfile(file):
        sys.exit(f'error: File "{file}" was not found.\n')
    return os.path.abspath(file)


def runSubprocess(cmd):
    """Run a subprocess

    Args:
        cmd (list): Command and args to run
    """
    sys.stderr.flush()
    cmd = list(map(str, cmd))
    print(f"Running {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    """if comp_proc.returncode != 0:
        sys.exit(f'''[{time.ctime()}] error: ProtHint exited due to an
                error in command: {" ".join(cmd)}''')"""


def systemCall(cmd):
    """Make a system call

    Args:
        cmd (string): Command to call
    """
    sys.stderr.flush()
    if subprocess.call(["bash", "-c", cmd]) != 0:
        sys.exit(f'''[{time.ctime()}] error: ProtHint exited due to an
                error in command: {cmd}''')


def callScript(name, args):
    """Call a script located in the ProtHint bin folder

    Args:
        name (string): Name of the script
        args (string): Command line arguments to use in the call
    """
    systemCall(binDir + '/' + name + ' ' + args)


def callDependency(name, args, location=''):
    """Call a dependency. ProtHint first looks for the dependency in the
    dependencies folder. If not found there, it tries to find it in the path.

    Args:
        name (string): Name of the dependency
        args (string): Command line arguments to use in the call
        location (string): Location of the dependency within the dependencies
                           folder
    """
    if location != '':
        location = location + '/'

    if os.path.isfile(binDir + '/../dependencies/' + location + name):
        systemCall(binDir + '/../dependencies/' + location + name + ' ' + args)
    else:
        sys.stderr.write('[' + time.ctime() + '] warning: Could not find ' +
                         name + ' in dependencies/' + location + ' folder.' +
                         ' Attempting to use ' + name + ' in the PATH.\n')
        if shutil.which(name):
            systemCall(f'{name} {args}')
        else:
            sys.exit(f'[{time.ctime()}] error: Could not find {name} in the PATH')


def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(
        description=f'''ProtHint {VERSION}: Pipeline for generating genome wide
         footprints of homologous proteins. The set of high confidence hints
         is generated using default thresholds in print_high_confidence.py
         script. If you wish to use different filtering criteria, re-run
         print_high_confidence.py script with custom thresholds.''')

    parser.add_argument('genome', metavar='genome.fasta', type=str,
                        help='Input genomic sequence in FASTA format.')
    parser.add_argument('proteins', metavar='proteins.fasta', type=str,
                        help='Protein database in FASTA format.')

    parser.add_argument('--workdir', type=str, default='.',
                        help='Folder for results and temporary files. If not specified, current directory is used')
    parser.add_argument('--geneSeeds', type=str, default='',
                        help='Gene prediction seeds in gtf format. If this file is provided, GeneMark-ES run is skipped.')
    parser.add_argument('--geneMarkGtf', type=str, default='',
                        help='Same as --geneSeeds, for backwards compatibility')
    parser.add_argument('--fungus', default=False, action='store_true',
                        help='Run GeneMark-ES in fungus mode.')
    parser.add_argument('--nonCanonicalSpliceSites', default=False, action='store_true',
                        help='Predict introns with non-canonical splices sites. By default, ProtHint only predicts introns with GT-AG, GC-AG, and AT-AC \
                        donor-acceptor sites.')
    parser.add_argument('--diamondPairs', type=str, default='',
                        help='File with "seed gene-protein" hits generated by DIAMOND. If this file is provided, DIAMOND search for protein hits is skipped.\
                        The seed genes in this file must correspond to seed genes passed by "--geneSeeds" option. All pairs in the file are used -- option \
                        "--maxProteinsPerSeed" is ignored.')
    parser.add_argument('--maxProteinsPerSeed', type=int, default=25,
                        help='Maximum number of protein hits per seed gene. Increasing this number leads to increased runtime and may improve the\
                        sensitivity of hints. Decreasing has an opposite effect. Default is set to 25.')
    parser.add_argument('--evalue', type=float, default=0.001,
                        help='Maximum e-value for DIAMOND alignments hits. Default = 0.001')
    parser.add_argument('--minExonScore', type=float, default=25,
                        help='Discard all hints inside/neighboring exons with score lower than minExonScore. Default = 25')
    parser.add_argument('--longGene', type=int, default=30000,
                        help='Threshold for what is considered a long gene in Spaln alignment. Default = 30000')
    parser.add_argument('--longProtein', type=int, default=15000,
                        help='Threshold for what is considered a long protein in Spaln alignment. Default = 15000')
    parser.add_argument('--cleanup', default=False, action='store_true',
                        help='Delete temporary files and intermediate results. Cleanup is turned off by default as it is useful to keep these files \
                        for troubleshooting and the intermediate results might be useful on their own.')
    parser.add_argument('--pbs', default=False, action='store_true',
                        help='Run GeneMark-ES and Spaln on pbs.')
    parser.add_argument('--threads', type=int, default=-1,
                        help='Number of threads used by ES, DIAMOND, and Spaln. By default, all available threads are used.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    parser.add_argument('--prevGeneSeeds', type=str, default='',
                        help='''File with gene seeds which were used in the
                        previous iteration. Next iteration of ProtHint is only
                        executed for --geneSeeds which differ from
                        --prevGeneSeeds. --prevSpalnGff is required when this
                        option is used since results from the previous
                        iteration are reused for seeds which do not differ
                        (Gene ids of such hints are updated to match the new
                         seed genes).''')
    parser.add_argument('--prevSpalnGff', type=str, default='',
                        help='Scored hints from previous iteration.')

    parser.add_argument('--ProSplign', default=False, action='store_true',
                        help=argparse.SUPPRESS)

    parser.add_argument('--diamond-db', default=False, type=str,
                        help='Prepared Diamond protein database.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
