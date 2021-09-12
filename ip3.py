import numpy as np
import random
import csv
import sys
import math
from itertools import product
from collections import defaultdict

# current set of tuning parameters and global variables is below
# Overlap hyperparam: haplotype with smaller Hamming distance will be chosen
overlap = 8
outfile = ""
end_bits = []

# return the transposed matrix
def reorient_genos(genos_t):
    gtlen = len(genos_t)
    slen = len(genos_t[0])
    genos = [[genos_t[jj][ii] for jj in range(gtlen)] for ii in range(slen)]
    return genos

# expects h_t to be the transpose of the matrix of haplotypes (i.e. h_t[0]
# is the first allele? of each haplotype
def print_haplotypes(h_t):
    with(open(outfile, "a+")) as f:
        for row in h_t:
            f.write(' '.join([str(r) for r in row]) + '\n')
    return

# expects a haplotype in the first of a list, i.e. ['1', '0', '1']
# returns the complementary haplotype, ['0', '1', '0']
def compl(haplotype):
    return [0 if h == 1 else 1 for h in haplotype]

# takes as input the genotype, a list, like ['1', '0', '0', '2', '1', '1'] and a
# permutation like ['0', '0', '1'], of length equal to the number of heterozygous
# alleles in the genotype and returns the corresponding haplotype in which
# 0's remain, 1's are mapped to 0's or 1's depending on the corresponding entry
# in the permutation and 2's are mapped to 1's
def gen_hapl(geno, perm, mperm):
    pidx = 0
    midx = 0
    glen = len(geno)
    hapl = []
    for ii in range(glen):
        if geno[ii] == 0:
            hapl.append(0)
        elif geno[ii] == 2:
            hapl.append(1)
        elif geno[ii] == -1:
            hapl.append(mperm[midx])
            midx += 1
        else:
            hapl.append(perm[pidx])
            pidx += 1
    return hapl

# for each genotype in genos, generates the set of all compatible haplotype
def gen_all_compatible(genos):
    haplotype_possibilities = []
    for geno in genos:
        haplo_temp = []
        hetero_cnt = geno.count(1)
        missing_cnt = geno.count(-1)
        if hetero_cnt == 0:
            if missing_cnt == 0:
                homo_hapl = tuple(gen_hapl(geno, [], []))
                haplotype_possibilities.append([(homo_hapl, homo_hapl)])
                continue
            for m in product([0, 1], repeat=missing_cnt):
                homo_hapl = tuple(gen_hapl(geno, [], m))
                haplo_temp.append((homo_hapl, homo_hapl))
            haplotype_possibilities.append(haplo_temp)
            continue
        for p in product([0, 1], repeat=hetero_cnt):
            # because we've been constructing the complement as well, by the time this
            # condition is satisfied, we will have examined all unique pairs of haplotypes
            if p[0] == 1:
                break

            if missing_cnt == 0:
                hapl = tuple(gen_hapl(geno, p, []))
                c_hapl = tuple(gen_hapl(geno, compl(p), []))
                haplo_temp.append((hapl, c_hapl))
                continue
            for m in product([0, 1], repeat=missing_cnt):
                hapl = tuple(gen_hapl(geno, p, m))
                c_hapl = tuple(gen_hapl(geno, compl(p), m))
                haplo_temp.append((hapl, c_hapl))
        haplotype_possibilities.append(haplo_temp)
    return haplotype_possibilities

# hamming distance between two haplotypes
def hd(hapl0, hapl1):
    return sum([h0 != h1 for h0, h1 in zip(hapl0, hapl1)])

# checks for parity between the end bits of each inferred haplotype and flips
# if necessary
def check_for_parity(hapl):
    num_hapls = len(hapl)
    for ii in range(0, num_hapls, 2):
        if hd(end_bits[ii], hapl[ii + 1][:overlap]) < hd(end_bits[ii], hapl[ii][:overlap]):
            temp_lst = hapl[ii + 1]
            hapl[ii + 1] = hapl[ii]
            hapl[ii] = temp_lst
    return

# updates the probability of each pair of haplotypes
def e_step(compat_probs, hapl_probs):
    for hapl_list in compat_probs:
        # two for loops; one just computes all P_{h_1}P_{h_2} and the other
        # updates the probabilities in compat_probs with P_{h_1}P_{h_2} / \sigma
        ph_sum = 0.0
        hlen = len(hapl_list)
        for ii in range(hlen):
            h1 = hapl_list[ii][0][0]
            h2 = hapl_list[ii][0][1]
            ph1_ph2 = hapl_probs[h1] * hapl_probs[h2]
            hapl_list[ii][1] = ph1_ph2
            ph_sum += ph1_ph2
        for ii in range(hlen):
            hapl_list[ii][1] /= ph_sum
    return

# updates the probability of each unique haplotype
def m_step(compat_probs, hapl_probs):
    two_n = float(2 * len(compat_probs))
    prob_dict = defaultdict(float)
    # sum all the entries in compat_probs corresponding to a given haplotype
    for hapl_list in compat_probs:
        for hapl_tuple in hapl_list:
            for hapl in hapl_tuple[0]:
                prob_dict[hapl] += hapl_tuple[1]
    # normalize the probabilities (divide by 2n) and assign to hapl_probs
    for uhapl in hapl_probs:
        hapl_probs[uhapl] = prob_dict[uhapl] / two_n
    return

def set_init_hapl_probs(hapl_set, hapl_probs):
    num_hapls = len(hapl_set)
    # initialize the probabilities for each unique haplotype
    for hapl in hapl_set:
        hapl_probs[hapl] = 1.0/num_hapls
    return

def q_metric(compat_probs, hapl_probs):
    q_sum = 0.0
    for hlist in compat_probs:
        max_prob = hlist[0][1]
        cur_pair = hlist[0][0]
        cnt_pairs = len(hlist)
        for ii in range(1, cnt_pairs):
            pair_prob = hlist[ii][1]
            if pair_prob > max_prob:
                max_prob = pair_prob
                cur_pair = hlist[ii][0]
        q_sum += np.log(hapl_probs[cur_pair[0]])
        q_sum += np.log(hapl_probs[cur_pair[1]])
    return q_sum

# function orchestrating the meat of the EM algorithm
def process_EM(genos):
    # for each genotype, generate a list of all pairs of compatible haplotypes
    compat_hapls = gen_all_compatible(genos)
    hapl_set = set()
    # generate the set of all unique haplotypes from the compatible haplotypes
    for hapl_list in compat_hapls:
        for hapl_tuple in hapl_list:
            for hapl in hapl_tuple:
                hapl_set.add(hapl)
    hapl_probs = {}
    set_init_hapl_probs(hapl_set, hapl_probs)
    compat_probs = []
    # initialize probabilities for each haplotype in C(g)
    for hapl_list in compat_hapls:
        aug_hapls = []
        l_len = float(len(hapl_list))
        for hapl_tuple in hapl_list:
            aug_hapls.append([hapl_tuple, 1/l_len])
        compat_probs.append(aug_hapls)
    # run the EM algorithm some number of times
    # could also fiddle with this hyperparameter
    for ii in range(6):
        m_step(compat_probs, hapl_probs)
        e_step(compat_probs, hapl_probs)
        #print(q_metric(compat_probs, hapl_probs))
    inferred_hapls = []
    for hlist in compat_probs:
        max_prob = hlist[0][1]
        cur_pair = hlist[0][0]
        cnt_pairs = len(hlist)
        for ii in range(1, cnt_pairs):
            pair_prob = hlist[ii][1]
            if pair_prob > max_prob:
                max_prob = pair_prob
                cur_pair = hlist[ii][0]
        inferred_hapls.append(cur_pair[0])
        inferred_hapls.append(cur_pair[1])
    global end_bits
    if len(end_bits) != 0:
        check_for_parity(inferred_hapls)
    h_t = reorient_genos(inferred_hapls)
    if len(end_bits) != 0:
        print_haplotypes(h_t[overlap:])
    else:
        print_haplotypes(h_t)
    end_bits = reorient_genos(h_t[-overlap:])
    return

def main(pfx):
    print(pfx)
    global outfile
    # Change masked_f depending on if you want to use an unmasked data file
    masked_f = pfx + '_masked.txt'
    #masked_f = pfx + '_gpred.txt'
    # change this line later to write to whatever the desired final output file is
    outfile = pfx + "_sol_wip3.txt"
    print(outfile)
    with open(outfile, 'w') as of:
        of.truncate()

    X = []
    conversions = {'0': 0, '1': 1, '2': 2, '*': -1}
    with(open(masked_f, 'r')) as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            X.append([int(conversions[ch]) for ch in row])
    num_snps = len(X)
    num_indivs = len(X[0])
    print('SNPs: {}'.format(num_snps))
    print('Individuals: {}'.format(num_indivs))

    # Fiddle with this hyperparam
    hetero_hyper = 14
    idx = 0
    num_hetero = [0 for x in range(num_indivs)]
    for ii in range(num_snps):
        for jj in range(num_indivs):
            if X[ii][jj] == 1 or X[ii][jj] == -1:
                num_hetero[jj] += 1
        if max(num_hetero) >= hetero_hyper:
            # No point in leaving the last couple SNPs lonely
            if ii < num_snps - 2:
                process_EM(reorient_genos(X[idx : ii]))
                idx = ii - overlap
                num_hetero = [0 for x in range(num_indivs)]
                # Need to fill num_hetero within the overlapped range
                for jj in range(idx, ii):
                    for kk in range(num_indivs):
                        if X[jj][kk] == 1 or X[jj][kk] == -1:
                            num_hetero[kk] += 1
                print(ii)
                continue
        if ii == num_snps - 1:
            process_EM(reorient_genos(X[idx:]))
            print(ii)

if __name__=="__main__":
    if len(sys.argv) != 2:
        print("Incorrect number of arguments. Provide one file prefix.")
    else:
        main(sys.argv[1])
