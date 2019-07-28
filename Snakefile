from __future__ import division, print_function
from string import ascii_uppercase
import itertools as itl
from math import log, ceil
from snakemake.utils import format as snakefmt
import random
import re
import json

configfile: "default_config.yml"

random.seed(int(config['seed']))

KMERLEN = config.get("kmerlen", 20)
N = config["population_size"]
COV_VAR = [(d['cov'], d['var']) for d in config["coverage_variablity"]]


#SEEDS = list(random.sample(range(10000), int(config["num_replicates"])))
SEEDS = set()
while len(SEEDS) < int(config["num_replicates"]):
    SEEDS.add(random.randrange(1e4))
SEEDS = list(SEEDS)
#print(SEEDS)

# Sample names
labels = ["".join(x)
          for x in itl.product(ascii_uppercase, repeat=int(ceil(log(N, 26))))]

GENOMES = [labels[i] for i in range(N)]
SAMPLES = [str(i+1) for i in range(config['runs_per_sample'])]


# All sim params
READ_NUMS = {}
SCALES = set()
for seed in SEEDS:
    random.seed(int(seed))
    readnum = {}
    for g in GENOMES:
        readnum[g] = {}
        for s in SAMPLES:
            readnum[g][str(s)] = {}
            for c, v in COV_VAR:
                cov = random.gauss(c, c * float(config["coverage_cv"]))
                nread = max((cov * config["genome_size"]) / 200, 100)
                try:
                    readnum[g][str(s)][str(c)][str(v)] = int(nread)
                except KeyError:
                    readnum[g][str(s)][str(c)] = {str(v): int(nread)}
                SCALES.add(v)
    READ_NUMS[str(seed)] = readnum
SCALES = list(SCALES)

METRICS = ['wip', 'ip']




rule all:
    input:
        expand(expand("data/{{seed}}/reads/{cov}x-{scale}/{{genome}}-{{sample}}_il.fastq.gz", zip,
                      cov=map(lambda x: x[0], COV_VAR),
                      scale=map(lambda x: x[1], COV_VAR)),
               seed=SEEDS, genome=GENOMES, sample=SAMPLES)


rule clean:
    shell:
        "rm -rf data .snakemake"


def scrm_args(wc):
    args = "{:d} {:d} ".format(config['population_size'], config['num_chroms'])
    #args += "-r {:d} {:d} ".format(Ne * Nrecomb, config['genome_size'])
    args += "-T "
    return args

rule population:
    output:
        "data/{seed}/population_prescale.nwk"
    params:
        args=scrm_args,
        seed=lambda w: w.seed,
    log:
        "data/{seed}/log/scrm.log"
    shell:
        "scrm"
        " {params.args}"
        " -seed {params.seed}"
        " 2>&1"
        "| tee {log}"
        "| grep '(' "  # Grep for (, which are only in a newick tree output line
        " >{output}"


rule scaletree:
    input:
        treefile="data/{seed}/population_prescale.nwk"
    output:
        treefile="data/{seed}/population.nwk"
    params:
        meandist=1.0,
    run:
        from ete3 import Tree
        import numpy as np
        from numpy import median, mean
        import itertools as itl
        from math import ceil, log
        from string import ascii_uppercase

        METRICS = [ 'mean', 'median', 'min', 'max', ]

        def pwdist(tree):
            '''Finds the (off-diagonal) pairwise distances between all tips of `tree`.
            '''
            dists = []
            for a, b in itl.combinations(tree.get_leaf_names(), 2):
                a = tree&a
                dists.append(a.get_distance(b))
            return np.array(dists)


        def normalise_tree(tree, to=1.0):
            '''
            Normalise branch lengths of `tree` such that mean of the pairwise
            distances is `to`.
            By default, normalise such that the mean of all pairwise distances is 1.0.
            '''
            dists = pwdist(tree)
            current = np.mean(dists)

            for node in tree.iter_descendants():
                node.dist /= current
            return tree


        def alphbetise_names(tree):
            '''Replace numeric tip labels with alphabetic ones. 1 -> A, 2 -> B etc.
            If there are more than 26 tips, labels are AA, AB, ..., ZZ and so forth for
            any number of tips.
            '''
            label_len = ceil(np.log(len(tree)) / np.log(26))  # how many letters do we need?
            labels = [''.join(letters)
                      for letters in itl.product(ascii_uppercase, repeat=label_len)]

            tiplabels = list(sorted(tree.get_leaf_names(), key=int))
            for i, leaf in enumerate(tiplabels):
                node = tree&leaf
                node.name = labels[i]
            return tree

        with open(input.treefile) as fh, open(output.treefile, "w") as ofh:
            for treeline in fh:
                tree = Tree(treeline)
                tree = alphbetise_names(tree)
                tree = normalise_tree(tree, params.meandist)
                print(tree.write(format=5), file=ofh)


rule seqgen:
    input:
        "data/{seed}/population.nwk"
    output:
        "data/{seed}/all_genomes-{scale}_seqgen.fasta"
    params:
        size=config["genome_size"],
    shell:
        "seq-gen"
        "   -m F84"
        "   -l {params.size}"
        "   -n 1"
        "   -f 0.3 0.2 0.2 0.3"
        "   -t 2.1"
        "   -of" # fasta
        "<{input} >{output}"


rule paramfile:
    input:
        "data/{seed}/population.nwk"
    output:
        "data/{seed}/dawg-{scale}.params",
    params:
        scale=lambda wc: wc.scale
    run:
        dawg = '''
        [Tree]
        Scale = {scale}
        [Indel]
        Model.Ins = Geo
        Model.Del = Geo
        Rate.Ins = 0.005
        Rate.Del = 0.005

        [Subst]
        Model = F84
        Freqs  = 0.3, 0.2, 0.2, 0.3
        Params = 2.5
        '''.format(scale=params.scale)
        with open(input[0])  as fh:
            for i, treeline in enumerate(fh):
                m = re.search(r'\[(.*)\](\(.*;)', treeline)
                if m:
                    length, tree = m.groups()
                else:
                    length = str(config['genome_size'])
                    tree = treeline
                dawg += '''
                [[-]]
                Root.Segment = {i:d}
                Root.Length = {length}
                Tree.Tree = {tree}
                '''.format(i=i, length=length, tree=tree)
        dawg = '\n'.join(x.lstrip() for x in dawg.split('\n'))
        with open(output[0], 'w') as fh:
            print(dawg, file=fh)


rule dawg:
    input:
        "data/{seed}/dawg-{scale}.params",
    output:
        "data/{seed}/all_genomes-{scale}_dawg.fasta"
    log:
        "data/{seed}/log/seqgen-{scale}.log"
    params:
        seed=lambda w: w.seed,
    shell:
        'dawg2'
        ' -o {output}'
        ' --seed {params.seed}'
        ' {input}'
        ' 2>{log} 1>&2'


rule alndist:
    input:
        ("data/{seed}/all_genomes-{scale}_dawg.fasta" if config["use_dawg"] else "data/{seed}/all_genomes-{scale}_seqgen.fasta")
    output:
        "data/{seed}/all_genomes-{scale}.dist"
    log:
        "data/{seed}/log/alndist-{scale}.log"
    run:
        from skbio import TreeNode, DistanceMatrix, TabularMSA, DNA
        from scipy.spatial.distance import hamming

        aln = TabularMSA.read(input[0], constructor=DNA)
        aln.reassign_index(minter="id")
        distmat = DistanceMatrix.from_iterable([seq.values for seq in aln],
                                            metric=hamming, keys=aln.index)
        distmat.write(output[0])


rule genomes:
    input:
        ("data/{seed}/all_genomes-{scale}_dawg.fasta" if config["use_dawg"] else 
        "data/{seed}/all_genomes-{scale}_seqgen.fasta")
    output:
        expand("data/{{seed}}/genomes/{{scale}}/{g}.fasta", g=GENOMES)
    params:
        dir="data/{seed}/genomes/{scale}/",
    log:
        "data/{seed}/log/splitfa-{scale}.log"
    run:
        import screed

        def seq2fa(name, seq, linelen=80):
            lines = ['>{}'.format(name), ]
            for start in range(0, len(seq), linelen):
                lines.append(seq[start:start + linelen])
            return '\n'.join(lines) + '\n'

        with screed.open(input[0]) as fh:
            for record in fh:
                name = record.name
                fname = "{}{}.fasta".format(params.dir, name)
                with open(fname, 'w') as ofh:
                    seq = str(record.sequence).translate({'-': ''})
                    print(seq2fa(name, seq, len(seq)), file=ofh)


rule mason_sim:
    input:
        "data/{seed}/genomes/{scale}/{genome}.fasta",
    output:
        r1=temp("data/{seed}/temp/{cov}x-{scale}/{genome}-{sample}_R1_mason.fastq"),
        r2=temp("data/{seed}/temp/{cov}x-{scale}/{genome}-{sample}_R2_mason.fastq"),
    params:
        rn=lambda w: str(READ_NUMS[w.seed][w.genome][w.sample][w.cov][w.scale]),
        seed=lambda w: w.seed,
    log:
        "data/{seed}/log/samples/{cov}x-{scale}/{genome}-{sample}.log"
    threads:
        1
    shell:
        "mason_simulator"
        " -ir {input}"
        " --illumina-read-length 101"
        " -o {output.r1}"
        " -or {output.r2}"
        " --seed {params.seed}"
        " -n {params.rn}"
        " --num-threads {threads}"
        " >{log} 2>&1"


rule wgsim:
    input:
        "data/{seed}/genomes/{scale}/{genome}.fasta",
    output:
        r1=temp("data/{seed}/temp/{cov}x-{scale}/{genome}-{sample}_R1_wgsim.fastq"),
        r2=temp("data/{seed}/temp/{cov}x-{scale}/{genome}-{sample}_R2_wgsim.fastq"),
    params:
        rn=lambda w: str(READ_NUMS[w.seed][w.genome][w.sample][w.cov][w.scale]),
        seed=lambda w: w.seed,
    log:
        "data/{seed}/log/samples/{cov}x-{scale}/{genome}-{sample}.log"
    threads:
        1
    shell:
        "wgsim"
        " -1 100" #r1 len
        " -2 100" #r2 len
        " -N {params.rn}"
        " -S {params.seed}"
        " {input}"
        " {output.r1}"
        " {output.r2}"
        " >{log} 2>&1"


rule ilfq:
    input:
        r1=("data/{seed}/temp/{cov}x-{scale}/{genome}-{sample}_R1_mason.fastq" if config.get("use_mason", True)
            else "data/{seed}/temp/{cov}x-{scale}/{genome}-{sample}_R1_wgsim.fastq"), 
        r2=("data/{seed}/temp/{cov}x-{scale}/{genome}-{sample}_R2_mason.fastq" if config.get("use_mason", True)
            else "data/{seed}/temp/{cov}x-{scale}/{genome}-{sample}_R2_wgsim.fastq"), 
    output:
        "data/{seed}/reads/{cov}x-{scale}/{genome}-{sample}_il.fastq.gz"
    log:
        "data/{seed}/log/join/{cov}x-{scale}/{genome}-{sample}.log"
    shell:
        "(seqhax pairs"
        "   {input.r1}"
        "   {input.r2}"
        "   -b >(gzip > {output})"
        ") 2>>{log}"

