import subprocess
import os
import sys
import time
import math
import random
from collections import Counter
import numpy as np

Max_iterations = 5
Wanted_ref_accuracy = 0.99
intToDNA = {'A': "A", 0: "A", 'C': "C", 1: "C", 'G': "G", 2: "G", 'T': "T", 3: "T"}


def Get_possible_refs(ref, cluster):
    n = len(ref)
    refs = [x for x in cluster if len(x) == n]
    refs.append(ref)
    return refs


"""ref is the reference strand"""
"""cluster is a list of strands"""
"""error rates has : substitution, insertion, deletion, similarity"""


def Align_with_ref(ref, cluster, error_rate, aligned_file_name):
    total = (error_rate[0] + error_rate[1] + error_rate[2]) * 100
    if total != 0:
        mis = -(math.log(100 * error_rate[0] + 1)) - 1
        ins = -(math.log(100 * error_rate[1] + 1)) - 1
        dele = -(math.log(100 * error_rate[2] + 1)) - 1
    else:
        mis = -1
        ins = -1
        dele = -1

    aligned_cluster = open(aligned_file_name, "w")
    max_len = len(max(cluster, key=len))
    deletion_hist = np.zeros((4, max_len))
    insertion_hist = np.zeros(max_len)
    substitution_hist = np.zeros(max_len)

    for strand in cluster:

        fasta1 = open("fasta1.txt", "w")
        fasta2 = open("fasta2.txt", "w")
        fasta1.write(">seq1\n")
        fasta1.write(ref)
        fasta2.write(">seq2\n")
        fasta2.write(strand)
        fasta1.close()
        fasta2.close()

        args=["./fogsaa.exe","fasta1.txt", "fasta2.txt", "1", "0", "0", str(mis), str(ins), str(dele), str(error_rate[3])]
        subprocess.call(args, stdin=None, stdout=None, stderr=None, shell=False)

        res_file1 = open("new_alseq1.txt", "r")
        res_file2 = open("new_alseq2.txt", "r")
        line1_aligned = res_file1.readline()
        line2_aligned = res_file2.readline()

        aligned_strand = ""

        """building an aligned strand with constant length"""
        for index, letter in enumerate(line1_aligned):
            if index >= max_len:
                deletion_hist=np.append(deletion_hist,[[0],[0],[0],[0]],axis=1)
                insertion_hist=np.append(insertion_hist,[0], axis=0)
                substitution_hist=np.append(substitution_hist,[0],axis=0)

            '''case of suspected deletion in the reference, in this case we delete the letter from the second strand'''
            if letter == "-":
                if line2_aligned[index] == 'A':
                    deletion_hist[0, index] += 1
                if line2_aligned[index] == 'C':
                    deletion_hist[1, index] += 1
                if line2_aligned[index] == 'G':
                    deletion_hist[2, index] += 1
                if line2_aligned[index] == 'T':
                    deletion_hist[3, index] += 1

                '''case of suspected insertion in the reference'''
            elif line2_aligned[index] == "-":
                aligned_strand += "-"
                insertion_hist[index] += 1

                '''case of substitution'''
            elif line2_aligned[index] != letter:
                aligned_strand += line2_aligned[index]
                substitution_hist[index] += 1

            else:
                aligned_strand += line2_aligned[index]

        aligned_cluster.write(aligned_strand)
        aligned_cluster.write("\n")

    fasta1.close()
    aligned_cluster.close()
    fasta2.close()
    res_file1.close()

    count = 0
    for i in range(len(ref)):
        if deletion_hist[0, i] >= len(cluster) / 2 or deletion_hist[1, i] >= len(cluster) / 2 or \
                deletion_hist[2, i] >= len(cluster) / 2 or deletion_hist[3, i] >= len(cluster) / 2 or \
                insertion_hist[i] >= len(cluster) / 2 or substitution_hist[i] >= len(cluster) / 2:
            count += 1
    evaluation = 1 - (count / len(ref))

    return evaluation, insertion_hist, deletion_hist, substitution_hist


def Align_without_ref(ref, cluster, error_rate, aligned_file_name):
    refs = Get_possible_refs(ref, cluster)

    random_ref = random.choice(refs)
    max_len = len(max(cluster, key=len))
    accuracy = 0
    count = 0

    used_ref=0
    used_insertion_hist = []
    used_deletion_hist = []
    used_substitution_hist = []

    while count <= Max_iterations:
        evaluation, insertion_hist, deletion_hist, substitution_hist = Align_with_ref(random_ref, cluster, error_rate,
                                                                                      aligned_file_name)
        if evaluation >= accuracy:
            used_ref = random_ref
            accuracy = evaluation
            used_insertion_hist = insertion_hist
            used_deletion_hist = deletion_hist
            used_substitution_hist = substitution_hist


        count += 1
        random_ref = random.choice(refs)



    return used_ref, accuracy, used_insertion_hist, used_deletion_hist, used_substitution_hist


def correct_ref(ref, aligned_cluster, insertion_hist, deletion_hist, substitution_hist):
    new_ref = ""
    np_aligned_cluster = np.array(aligned_cluster)
    ref_index = 0
    hist_index = 0
    while hist_index < len(substitution_hist):

        if substitution_hist[hist_index] > (len(aligned_cluster) / 2):
            options = np.array([strand[ref_index] for strand in np_aligned_cluster])
            c = Counter(options)
            value, num = c.most_common()[0]
            if value == "-":
                continue
            new_ref += value
            ref_index += 1

        elif deletion_hist[0, hist_index] > (len(np_aligned_cluster) / 2):
            new_ref += "A"
        elif deletion_hist[1, hist_index] > (len(np_aligned_cluster) / 2):
            new_ref += "C"
        elif deletion_hist[2, hist_index] > (len(np_aligned_cluster) / 2):
            new_ref += "G"
        elif deletion_hist[3, hist_index] > (len(np_aligned_cluster) / 2):
            new_ref += "T"

        elif insertion_hist[hist_index] > (len(aligned_cluster) / 2):
            ref_index += 1

        elif ref_index < len(ref):
            new_ref += ref[ref_index]
            ref_index += 1
        hist_index += 1

    return new_ref


def Align(ref, cluster, error_rate, with_ref, enable_two_rounds):
    return_file = "firstPhaseAlignment.txt"
    used_ref = ""
    evaluation = 0
    deletion_hist = []

    if with_ref == True:
        used_ref = ref
        evaluation, insertion_hist, deletion_hist, substitution_hist = Align_with_ref(ref, cluster, error_rate,
                                                                                      "firstPhaseAlignment.txt")
    else:
        used_ref, evaluation, insertion_hist, deletion_hist, substitution_hist = Align_without_ref(ref, cluster,
                                                                                                   error_rate,
                                                                                                   "firstPhaseAlignment.txt")

    if enable_two_rounds == True:
        if evaluation < Wanted_ref_accuracy:
            new_ref = correct_ref(used_ref, cluster, insertion_hist, deletion_hist, substitution_hist)
            new_evaluate, new_insertion_hist, new_deletion_hist, new_substitution_hist = Align_with_ref(new_ref,
                                                                                                        cluster,
                                                                                                        error_rate,
                                                                                                        "SecondPhaseAlignment.txt")

            if new_evaluate > evaluation:
                used_ref = new_ref
                evaluation = new_evaluate
                insertion_hist = new_insertion_hist
                deletion_hist = new_deletion_hist
                substitution_hist = new_substitution_hist
                return_file = "SecondPhaseAlignment.txt"

    return used_ref, evaluation, insertion_hist, deletion_hist, substitution_hist, return_file


def get_pred(pred_probs, length):
    pred = ""
    for k in range(length):
        pred_index = [pred_probs[j][k] for j in [0, 1, 2, 3]]
        pred_val = np.argmax(pred_index)  # , key=map_index)
        pred = pred + intToDNA[pred_val]
    return pred

#command line: python main.py cluster_file_path substitution_rate insertion_rate deletion_rate similarity with_ref 2_round
#similarity = 1-subs-ins-dele
#with_ref is a boolean, to allow using the reference or not
#2_round is a boolean, to allow the improvement round
if __name__ == '__main__':
    args = sys.argv
    if len(args) != 8:
        print("Invalid Arguments")

    else:
        input = open(args[1])
        lines = input.readlines()
        ref = lines[0]
        cluster = lines[1:]

        error_rate = [float(args[2]), float(args[3]), float(args[4]), float(args[5])]
        if args[6] == "True":
            with_ref = True
        else:
            with_ref = False
        if args[7] == "True":
            second_round = True
        else:
            second_round = False

        used_ref, evaluation, insertion_hist, deletion_hist, substitution_hist, return_file = Align(ref, cluster, error_rate, with_ref, second_round)

        print("output file name: " + return_file)



    """testing the code"""
    """ref = "CACGCTCAGCATGCGACATGCTCGCGCAGCACGTGTAGCTGTGTCTCACTGCTCTGCTACGAGTGCTAGTGCACACGTGCGCATCTACACAGATATAGACGTGAGCTGTAGTGCTACAGTCGTGCTCAGCAGATGTCTACACGCAGTAGTCATGACACGCTAGTGCACGAGCACGCGATGTCGTCGACTGTAGACGCACACGCGCGTATAGAGACATAGC"
    cluster = ["TCACCAGGTCAAATCCCACCAATAACATGCGCAGACGTGGCGACGATACCCTTGATTGGTGCTCCACCTTGTCTCCGTTGTCTATCGAGTTGCTTAAGTTCGCCATCGCGCGTCATTCTTACAGCAAGATATAGAATCGTCGCCAGTATTGGCTACAAGTCGTAGCTTCAAGTTAGGACTGTTCTACACAGCAGTTCGTCCAGGCTACCGCTCGTGCAACGAGCCAACGCGATGTTCCGTCGACGTGTACGCCGACGCTGCTGTTAGAAAGTATTAGC"
               ,"TCCAGGCATCAGTCCTACTGCAGCGGTGTCGGCTACGTGTCTTTTCTTGCACCTAGGGTGCCTTAACGTGCAGGCAGACTGTAGCGAAGTCTACTCAGATGATAGAACGGGATGCGGTAGTGCGACAGTCCGTACCTCAGATGATTGTGCTACTATCAGCATGATTAGTGCATCGAACCACGTCTTGTTGACAGCGAGCACAGCCGAGGTGCGTCGGTCTGCTGAGAACCGGCACACCGGTGTATTAAGACGACATAACC",
               "ACCATAACCTTCTACCTACTGACGTAAAATTTGACTACAGTCTATAATCAGCGAGAGGTCCTGCCCTCCTACCATCCCGGCTTCTACAGAAAGTGCTCGTGCACACGTGACTGTCATCTAGCACAGATGATACATGTAGAGCTGTACGTGGTTCCAGTCCTGTGCTACGAGCGATGTCTGCGAACGGCAAGTAGGTGTACTGACATCGACTTCCAACGAGCCACAGCGTATGTCGCTGCGTCTGTGAGCACGCCCACGCCGTCGACGTAGTGAGGAGGGGGATTAGGC",
               "TCTCTTCTCCAACCGTCCTGCGACTATCGGTAACCTCGGATGATAGTCAGGAGATCCTTCAACTTAGTCATATTGCGCTGAATGTTGGCTAGGCTTGGCGAGCAATCGATGGCGGCCATACTGACAACAAGAATAATTAGAGGTGGCTAGTTGAGTTGGCTACAGTCGTGCCTTCGAAACATCGGTCCTGACACGCTAGTCAGTAACCTACTAGTGTCCGAGCAACGACGATGTTCGTCGACCTGTAGCCGCATAAGACGGCGTTCACTGAGAGAGCATCGGC",
               "CACGCTCAGCATGCGACATGGTCGCGCAGCACGTGTAGCTGTGTCTCACTGCTCTGCTACGTTTGCTAGTGCACACGTGCGCATCTACACAGATTTAGACGTGGGCTGTAGTGCTACAGTCGTGCTCAGCAGATGTCTACACGCAGTAGTCATGACACGCTAGTGCACGAGCACGCGATGTCGTCGACTGTAGACGCACACGCGCGTATAGAGACATAGC"]
    error_rate = [0.121,0.367,0.0433,0.47]

    used_ref, evaluation, insertion_hist, deletion_hist, substitution_hist, return_file = Align(ref,cluster,error_rate,False,True)

    print("evaluate: " + str(evaluation))"""









    """# Path to per-cluster results
    # cluster_results_path = "./../../dgx_data/results/pilot_nanopore_adaptor_10_r/pseudo"
    cluster_results_path = "/Users/omersabary/Downloads/15.05.2022.pilot_nanopore_full_test_1/pseudo/pilot_nanopore.DNAformer_s_v1.ce_sl1.size_32/pilot_nanopore_full_test_1/cluster_size_30/"
    results = {"test_1": {}, "test_2": {}, "test_3": {}}
    folder_names = [f for f in listdir(cluster_results_path)]
    print(folder_names)
    folder_names.sort()
    read_length = 128
    num_of_clusters = 0
    for f_name in ["test_1"]:  # folder_names[:1]:
        if not f_name.startswith("t"): continue
        # name = f_name.split('_')
        # size = 'size_' + name[-1]
        folder_path = cluster_results_path + '/' + f_name + '/'
        test_name = f_name

        pred_files = [f for f in listdir(folder_path)]
        pred_files.sort()
        num_of_clusters = len(pred_files)
        pred_files_short = pred_files
        for json_file in tqdm(pred_files_short):
            path = folder_path + '/' + json_file
            with open(path) as f:
                data = json.load(f)
                cluster_size = len(data["noisy_copies"])
                label = getDNAstrand(data["label"])
                pred = get_pred(pred_probs=data["pred_probs"], length=len(label))
                res = (label == pred)
                dist = hamming_distance(pred, label)
                pred_left = get_pred(pred_probs=data["pred_left"], length=len(label))
                pred_right = get_pred(pred_probs=data["pred_right"], length=len(label))
                results[test_name][data["index"]] = {'label': label, 'pred': pred, 'cluster_size': cluster_size,
                                                     'dist': dist, 'res': res, 'pred_left': pred_left,
                                                     'pred_right': pred_right, 'pred_probs': data["pred_probs"]}"""


