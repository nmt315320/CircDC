__author__ = 'Fule Liu'

import sys
import os
import pickle
from math import pow
import time

import const
from util import frequency
from util import get_data
from util import check_args, read_k, write_to_file
from nac import make_kmer_list
from data import index_list


"""Prepare for PseKNC."""


class AAIndex():
    def __init__(self, head, index_dict):
        self.head = head
        self.index_dict = index_dict

    def __str__(self):
        return "%s\n%s" % (self.head, self.index_dict)


def pseknc(input_data, k, w, lamada, phyche_list, alphabet, extra_index_file=None, all_prop=False, theta_type=1):
    """This is a complete process in PseKNC.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param extra_index_file: a file path includes the user-defined phyche_index.
    :param all_prop: bool, choose all physicochemical properties or not.
    """
    phyche_list = get_phyche_list(k, phyche_list,
                                  extra_index_file=extra_index_file, alphabet=alphabet, all_prop=all_prop)
    # Get phyche_vals.
    if alphabet == index_list.DNA or alphabet == index_list.RNA:
        if extra_index_file is not None:
            extra_phyche_index = get_extra_index(extra_index_file)
            from util import normalize_index

            phyche_vals = get_phyche_value(k, phyche_list, alphabet,
                                           normalize_index(extra_phyche_index, alphabet, is_convert_dict=True))
        else:
            phyche_vals = get_phyche_value(k, phyche_list, alphabet)
    elif alphabet == index_list.PROTEIN:
        phyche_vals = get_aaindex(phyche_list)
        if extra_index_file is not None:
            phyche_vals.extend(extend_aaindex(extra_index_file))

    seq_list = get_data(input_data, alphabet)

    return make_pseknc_vector(seq_list, phyche_vals, k, w, lamada, alphabet, theta_type)


def ipseknc(input_data, k, w, lamada, phyche_list, alphabet, extra_index_file=None, all_prop=False):
    """This is a complete process in iPseKNC, k is kmer, but the index is just for dinucleotide.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param extra_index_file: a file path includes the user-defined phyche_index.
    :param all_prop: bool, choose all physicochemical properties or not.
    """
    phyche_list = get_phyche_list(k=2, phyche_list=phyche_list,
                                  extra_index_file=extra_index_file, alphabet=alphabet, all_prop=all_prop)

    # Get phyche_vals.
    if extra_index_file is not None:
        extra_phyche_index = get_extra_index(extra_index_file)
        from util import normalize_index

        phyche_vals = get_phyche_value(k=2, phyche_list=phyche_list, alphabet=alphabet,
                                       extra_phyche_index=normalize_index(extra_phyche_index, alphabet,
                                                                          is_convert_dict=True))
    else:
        phyche_vals = get_phyche_value(k=2, phyche_list=phyche_list, alphabet=alphabet)

    seq_list = get_data(input_data, alphabet)

    return make_pseknc_vector(seq_list, phyche_vals, k, w, lamada, alphabet, theta_type=3)


def get_phyche_list(k, phyche_list, extra_index_file, alphabet, all_prop=False):
    """Get phyche_list and check it.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param all_prop: bool, choose all physicochemical properties or not.
    """
    if phyche_list is None or len(phyche_list) == 0:
        if extra_index_file is None and all_prop is False:
            error_info = 'Error, The phyche_list, extra_index_file and all_prop can\'t be all False.'
            raise ValueError(error_info)

    from data import index_list

    # Set all_prop_list.
    all_prop_list = []
    try:
        if alphabet == index_list.DNA:
            if k == 2:
                all_prop_list = index_list.didna_list
            elif k == 3:
                all_prop_list = index_list.tridna_list
            else:
                error_info = 'Error, the k value must be 2 or 3.'
                raise ValueError(error_info)
        elif alphabet == index_list.RNA:
            if k == 2:
                all_prop_list = index_list.dirna_list
            else:
                error_info = 'Error, the k or alphabet error.'
                raise ValueError(error_info)
        elif alphabet == index_list.PROTEIN:
            all_prop_list = index_list.pro_list
        else:
            error_info = "Error, the alphabet must be dna, rna or protein."
            raise ValueError(error_info)
    except:
        raise

    # Set and check physicochemical properties.
    try:
        # Set all properties.
        if all_prop is True:
            phyche_list = all_prop_list
        # Check phyche properties.
        else:
            for e in phyche_list:
                if e not in all_prop_list:
                    error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                    raise NameError(error_info)
    except:
        raise

    return phyche_list


def get_extra_index(filename):
    """Get the extend indices from index file, only work for DNA and RNA."""
    extra_index_vals = []
    with open(filename) as f:
        lines = f.readlines()
        for ind, line in enumerate(lines):
            if line[0] == '>':
                vals = lines[ind + 2].rstrip().split('\t')
                vals = [float(val) for val in vals]
                extra_index_vals.append(vals)

    return extra_index_vals


def get_aaindex(index_list):
    """Get the aaindex from data/aaindex.data.

    :param index_list: the index we want to get.
    :return: a list of AAIndex obj.
    """
    new_aaindex = []
    with open('data/aaindex.data', 'rb') as f:
        aaindex = pickle.load(f)
        for index_vals in aaindex:
            if index_vals.head in index_list:
                new_aaindex.append(index_vals)

    return new_aaindex


def extend_aaindex(filename):
    """Extend the user-defined AAIndex from user's file.
    :return: a list of AAIndex obj.
    """
    from scrip.extract_aaindex import extra_aaindex, norm_index_vals

    aaindex = extra_aaindex(filename)
    for ind, e in enumerate(aaindex):
        aaindex[ind] = AAIndex(e.head, norm_index_vals(e.index_dict))

    return aaindex


def get_ext_ind_pro(filename):
    """Get the extend indices from index file, only work for protein."""
    inds = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    aaindex = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line[0] == '>':
                temp_name = line[1:].rstrip()
                vals = lines[i + 2].rstrip().split('\t')
                ind_val = {ind: float(val) for ind, val in zip(inds, vals)}
                aaindex.append(AAIndex(temp_name, ind_val))

    return aaindex


def get_phyche_value(k, phyche_list, alphabet, extra_phyche_index=None):
    """Generate DNA or RNA phyche_value.

    :param k: int, the value of k-tuple.
    :param phyche_list: physicochemical properties list.
    :param extra_phyche_index: dict, the key is the olinucleotide (string),
                                     the value is its physicochemical property value (list).
                               It means the user-defined physicochemical indices.
    """
    if extra_phyche_index is None:
        extra_phyche_index = {}

    phyche_value = extend_phyche_index(get_phyche_index(k, phyche_list, alphabet), extra_phyche_index)

    return phyche_value


def extend_phyche_index(original_index, extend_index):
    """Extend DNA or RNA {phyche:[value, ... ]}"""
    if extend_index is None or len(extend_index) == 0:
        return original_index
    for key in list(original_index.keys()):
        original_index[key].extend(extend_index[key])
    return original_index


def get_phyche_factor_dic(k, alphabet):
    """Get all DNA or RNA {nucleotide: [(phyche, value), ...]} dict."""
    full_path = os.path.realpath(__file__)
    if 2 == k and alphabet == index_list.DNA:
        file_path = "%s/data/didna.data" % os.path.dirname(full_path)
    elif 2 == k and alphabet == index_list.RNA:
        file_path = "%s/data/dirna.data" % os.path.dirname(full_path)
    elif 3 == k:
        file_path = "%s/data/mmc4.data" % os.path.dirname(full_path)
    else:
        sys.stderr.write("The k can just be 2 or 3.")
        sys.exit(0)

    try:
        with open(file_path, 'rb') as f:
            phyche_factor_dic = pickle.load(f)
    except:
        with open(file_path, 'r') as f:
            phyche_factor_dic = pickle.load(f)

    return phyche_factor_dic


def get_phyche_index(k, phyche_list, alphabet):
    """get phyche_value according phyche_list."""
    phyche_value = {}
    if 0 == len(phyche_list):
        for nucleotide in make_kmer_list(k, alphabet):
            phyche_value[nucleotide] = []
        return phyche_value

    nucleotide_phyche_value = get_phyche_factor_dic(k, alphabet)
    for nucleotide in make_kmer_list(k, alphabet):
        if nucleotide not in phyche_value:
            phyche_value[nucleotide] = []
        for e in nucleotide_phyche_value[nucleotide]:
            if e[0] in phyche_list:
                phyche_value[nucleotide].append(e[1])

    return phyche_value


"""Calculate PseKNC."""


def parallel_cor_function(nucleotide1, nucleotide2, phyche_index):
    """Get the cFactor.(Type1)"""
    temp_sum = 0.0
    phyche_index_values = list(phyche_index.values())
    len_phyche_index = len(phyche_index_values[0])
    for u in range(len_phyche_index):
        temp_sum += pow(float(phyche_index[nucleotide1][u]) - float(phyche_index[nucleotide2][u]), 2)

    return temp_sum / len_phyche_index


def series_cor_function(nucleotide1, nucleotide2, big_lamada, phyche_value):
    """Get the series correlation Factor(Type 2)."""
    return float(phyche_value[nucleotide1][big_lamada]) * float(phyche_value[nucleotide2][big_lamada])


def pro_cor_fun1(ri, rj, aaindex_list):
    _sum = 0.0
    len_index = len(aaindex_list)
    for aaindex in aaindex_list:
        _sum += pow(aaindex.index_dict[ri] - aaindex.index_dict[rj], 2)
    return _sum / len_index


def pro_cor_fun2(ri, rj, aaindex):
    return aaindex.index_dict[ri] * aaindex.index_dict[rj]


def get_parallel_factor(k, lamada, sequence, phyche_value, alphabet):
    """Get the corresponding factor theta list."""
    theta = []
    l = len(sequence)

    for i in range(1, lamada + 1):
        temp_sum = 0.0
        for j in range(0, l - k - i + 1):
            nucleotide1 = sequence[j: j + k]
            nucleotide2 = sequence[j + i: j + i + k]
            if alphabet == index_list.DNA or alphabet == index_list.RNA:
                temp_sum += parallel_cor_function(nucleotide1, nucleotide2, phyche_value)
            elif alphabet == index_list.PROTEIN:
                temp_sum += pro_cor_fun1(nucleotide1, nucleotide2, phyche_value)

        theta.append(temp_sum / (l - k - i + 1))

    return theta


def get_series_factor(k, lamada, sequence, phyche_value, alphabet):
    """Get the corresponding series factor theta list."""
    theta = []
    l_seq = len(sequence)
    if alphabet == index_list.DNA or alphabet == index_list.RNA:
        temp_values = list(phyche_value.values())
        max_big_lamada = len(temp_values[0])
    elif alphabet == index_list.PROTEIN:
        max_big_lamada = len(phyche_value)

    for small_lamada in range(1, lamada + 1):
        for big_lamada in range(max_big_lamada):
            temp_sum = 0.0
            for i in range(0, l_seq - k - small_lamada + 1):
                nucleotide1 = sequence[i: i + k]
                nucleotide2 = sequence[i + small_lamada: i + small_lamada + k]
                if alphabet == index_list.DNA or alphabet == index_list.RNA:
                    temp_sum += series_cor_function(nucleotide1, nucleotide2, big_lamada, phyche_value)
                elif alphabet == index_list.PROTEIN:
                    temp_sum += pro_cor_fun2(nucleotide1, nucleotide2, phyche_value[big_lamada])

            theta.append(temp_sum / (l_seq - k - small_lamada + 1))

    return theta


def make_pseknc_vector(sequence_list, phyche_value, k=2, w=0.05, lamada=1, alphabet=index_list.DNA, theta_type=1):
    """Generate the pseknc vector."""
    kmer = make_kmer_list(k, alphabet)
    vector = []

    for sequence in sequence_list:
        if len(sequence) < k or lamada + k > len(sequence):
            error_info = "Sorry, the sequence length must be larger than " + str(lamada + k)
            sys.stderr.write(error_info)
            sys.exit(0)

        # Get the nucleotide frequency in the DNA sequence.
        fre_list = [frequency(sequence, str(key)) for key in kmer]
        fre_sum = float(sum(fre_list))

        # Get the normalized occurrence frequency of nucleotide in the DNA sequence.
        fre_list = [e / fre_sum for e in fre_list]

        # Get the theta_list.
        if 1 == theta_type:
            theta_list = get_parallel_factor(k, lamada, sequence, phyche_value, alphabet)
        elif 2 == theta_type:
            theta_list = get_series_factor(k, lamada, sequence, phyche_value, alphabet)
        elif 3 == theta_type:
            theta_list = get_parallel_factor(k=2, lamada=lamada, sequence=sequence,
                                             phyche_value=phyche_value, alphabet=alphabet)
        theta_sum = sum(theta_list)

        # Generate the vector according the Equation 9.
        denominator = 1 + w * theta_sum

        temp_vec = [round(f / denominator, 8) for f in fre_list]
        for theta in theta_list:
            temp_vec.append(round(w * theta / denominator, 8))

        vector.append(temp_vec)

    return vector


def read_index(index_file):
    with open(index_file) as f_ind:
        lines = f_ind.readlines()
        ind_list = [index.rstrip() for index in lines]
        return ind_list


def main(args):

    file_list = args.inputfiles
    outputfile_list = []
    label_list = args.labels
    output_format = args.f
    if len(file_list) == 0:
        print 'Input files not found.'
        return False
    if output_format == 'svm' and len(label_list) == 0:
        print 'The labels of the input files should be set.'
        return False
    if output_format == 'svm' and len(file_list) != len(label_list):
        print 'The number of labels should be the same as that of the input files.'
        return False

    if args.out is not None:
        outputfile_list = args.out
        if len(outputfile_list) != len(file_list):
            print 'The number of output files should be the same as that of input files.'
            return False
    elif args.out is None:
        outputfile_list = []
        if output_format =='svm':
            for in_file_name in file_list:
                file_elem_list = list(os.path.splitext(in_file_name))
                out_name = file_elem_list[0] + '_svm' + file_elem_list[1]
                outputfile_list.append(out_name)
        elif output_format =='tab':
            for in_file_name in file_list:
                file_elem_list = list(os.path.splitext(in_file_name))
                out_name = file_elem_list[0] + '_tab' + file_elem_list[1]
                outputfile_list.append(out_name)
        elif output_format =='csv':
            for in_file_name in file_list:
                file_elem_list = list(os.path.splitext(in_file_name))
                out_name = file_elem_list[0] + '_csv' + file_elem_list[1]
                outputfile_list.append(out_name)
    if output_format !='svm':
        label_list = [0] * len(file_list)

    for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
        with open(input_file) as f:
            # Get index_list.
            if args.i is not None:
                ind_list = read_index(args.i)
            else:
                ind_list = []

            default_e = []
            # Set Pse default index_list.
            if args.alphabet == 'DNA':
                alphabet_list = index_list.DNA
                if args.k == 2:
                    default_e = const.DI_INDS_6_DNA
                elif args.k == 3:
                    default_e = const.TRI_INDS_DNA
            elif args.alphabet == 'RNA':
                alphabet_list = index_list.RNA
                default_e = const.DI_INDS_RNA
            elif args.alphabet == 'Protein':
                alphabet_list = index_list.PROTEIN
                default_e = const.INDS_3_PROTEIN

            theta_type = 1
            if args.method in const.THETA_1_METHODS:
                theta_type = 1
            elif args.method in const.THETA_2_METHODS:
                theta_type = 2
            elif args.method == 'PseKNC':
                theta_type = 3
            else:
                print("Method error!")

            # PseKNC.
            if args.method != 'PseKNC':
                if args.e is None and len(ind_list) == 0 and args.a is False:
                    # Default Pse.
                    res = pseknc(f, args.k, args.w, args.lamada, default_e, alphabet_list,
                                 extra_index_file=args.e, all_prop=args.a, theta_type=theta_type)
                else:
                    res = pseknc(f, args.k, args.w, args.lamada, ind_list, alphabet_list,
                                 extra_index_file=args.e, all_prop=args.a, theta_type=theta_type)
            # iPseKNC.
            else:
                if args.e is None and len(ind_list) == 0 and args.a is False:
                    # Default iPse.
                    res = ipseknc(f, args.k, args.w, args.lamada, const.DI_INDS_6_DNA, alphabet_list,
                                  extra_index_file=args.e, all_prop=args.a)
                else:
                    res = ipseknc(f, args.k, args.w, args.lamada, ind_list, alphabet_list,
                                  extra_index_file=args.e, all_prop=args.a)
        write_to_file(res, output_format, label, output_file)

    if len(outputfile_list) != 0:
        for index, output_file in enumerate(outputfile_list):
            out_with_full_path = os.path.abspath(output_file)
            if os.path.isfile(out_with_full_path):
                if index == 0:
                    print 'The output file(s) can be found here:'
                print out_with_full_path



if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter

    parse = argparse.ArgumentParser(description="This is pse module for generate pse vector.",
                                    formatter_class=RawTextHelpFormatter)

    parse.add_argument('inputfiles', nargs='*',
                       help="The input files in FASTA format.")
    parse.add_argument('-out', nargs='*',
                       help="The output files for storing feature vectors.")
    parse.add_argument('alphabet', choices=['DNA', 'RNA', 'Protein'],
                       help="The alphabet of sequences.")
    parse.add_argument('method', type=str,
                       help="The method name of pseudo components.")

    parse.add_argument('-lamada', type=int, default=2,
                       help="The value of lamada. default=2")
    parse.add_argument('-w', type=float, default=0.1,
                       help="The value of weight. default=0.1")
    parse.add_argument('-k', type=int,
                       help="The value of kmer, it works only with PseKNC method.")
    parse.add_argument('-i',
                       help="The indices file user choose.\n"
                            "Default indices:\n"
                            "DNA dinucleotide: Rise, Roll, Shift, Slide, Tilt, Twist.\n"
                            "DNA trinucleotide: Dnase I, Bendability (DNAse).\n"
                            "RNA: Rise, Roll, Shift, Slide, Tilt, Twist.\n"
                            "Protein: Hydrophobicity, Hydrophilicity, Mass.")
    parse.add_argument('-e', help="The user-defined indices file.\n")
    parse.add_argument('-all_index', dest='a', action='store_true', help="Choose all physicochemical indices")
    parse.add_argument('-no_all_index', dest='a', action='store_false',
                       help="Do not choose all physicochemical indices, default.")
    parse.set_defaults(a=False)
    parse.add_argument('-f', default='tab', choices=['tab', 'svm', 'csv'],
                       help="The output format (default = tab).\n"
                            "tab -- Simple format, delimited by TAB.\n"
                            "svm -- The libSVM training data format.\n"
                            "csv -- The format that can be loaded into a spreadsheet program.")

    parse.add_argument('-labels', nargs='*',
                       help="The labels of the input files.\n"
                       "For binary classification problem, the labels can only be '+1' or '-1'.\n"
                       "For multiclass classification problem, the labels can be set as a list of integers.")


    args = parse.parse_args()
    args.k = read_k(args.alphabet, args.method, args.k)

    # print(args)
    if check_args(args, 'pse.py'):
        print("Calculating...")
    start_time = time.time()
    main(args)
    print("Done.")
    print("Used time: %.2fs" % (time.time() - start_time))
