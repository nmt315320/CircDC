__author__ = 'Fule Liu'

import time, re, os, sys, stat
import subprocess
import const
from util import get_data, check_args, read_k, write_to_file
from pse import get_phyche_list, get_extra_index, get_phyche_value, get_aaindex, extend_aaindex, AAIndex
from data import index_list


# ==============================================================================
# autocorrelation: moreau,geary,moran
def getValues(prop, supInfo):
    values = ""
    name = re.search(prop, supInfo)
    if name:
        strr = prop + '\s*\,(.+)'
        b = re.search(strr, supInfo)
        if b:
            values = b.group(1)
    return values


def sepSequence(seq, k):
    i = k - 1
    seqq = []
    while i < len(seq):
        j = 0
        nuc = ''
        while j < k:
            nuc = seq[i - j] + nuc
            j = j + 1
        seqq.append(nuc)
        i += 1
    return seqq


def getSpecificValue(olinuc, olinucs, prop, supInfo):
    # olinucs = getNucSeq(SupFileName).split(",")
    olinucs = olinucs.split(",")
    values = getValues(prop, supInfo).rstrip()
    values = values.split(",")
    # valueS = [float(x) for x in values.split(",")]
    count = olinucs.index(olinuc)
    value = values[count]
    return float(value)


def avgP(seq, olinucs, length, k, prop, supInfo):
    limit = length - k + 1
    i = 1
    sum = 0
    while i < limit or i == limit:
        # value = hn(seq[i - 1],prop,SupFileName)
        value = getSpecificValue(seq[i - 1], olinucs, prop, supInfo)
        sum = sum + value
        i = i + 1
    sum = sum / limit
    return sum


#  geary
# --------------------------------------
# inputs: seq = string, length = int, k = int, l = int, prop = string,
#         SupFileName = string
# output: final = int
def geary(seq, olinucs, length, k, l, prop, supInfo):
    lim = length - k + 1
    limit = length - k - l + 1
    b = 1
    sqr = 0
    while b < limit or b == limit:
        current = getSpecificValue(seq[b - 1], olinucs, prop, supInfo)
        # hn(seq[b-1],prop,SupFileName)
        next = getSpecificValue(seq[b + l - 1], olinucs, prop, supInfo)
        # hn(seq[b+l-1],prop,SupFileName)
        sqr = sqr + ((current - next) * (current - next))
        b = b + 1
    top = sqr * lim
    limit2 = (length - k - l + 1)
    c = 1
    sqr2 = 0
    while c < limit2 or c == limit2:
        # current = hn(seq[c-1],prop,SupFileName)
        current = getSpecificValue(seq[c - 1], olinucs, prop, supInfo)
        avg = avgP(seq, olinucs, length, k, prop, supInfo)
        sqr2 = sqr2 + (current - avg) * (current - avg)
        c = c + 1
    bottom = sqr2 * limit * 2
    final = float((top / bottom) * 1000) / 1000.0
    return final


#  Moreau
# -------------------------------------
# inputs: seq = string, length = int, k = int, l = int, prop = string,
#         supFileName = string
# output: final = int

def moreau(seq, olinucs, length, k, l, prop, supInfo):
    limit = length - k - l + 1
    d = 1
    prod = 0
    while d < limit or d == limit:
        current = getSpecificValue(seq[d - 1], olinucs, prop, supInfo)
        # hn(seq[d-1],prop,SupFileName)
        next = getSpecificValue(seq[d + l - 1], olinucs, prop, supInfo)
        # hn(seq[d+l-1],prop,SupFileName)
        prod = prod + (current * next)
        d = d + 1
    final = prod / limit
    return final

#  moran
# --------------------------------------
# inputs: seq = string, length = int, k = int, l = int, prop = string,
#         SupFileName = string
# output: final = int

def moran(seq, olinucs, length, k, l, prop, supInfo):
    limit = length - k - l + 1
    j = 1
    top = 0
    avg = avgP(seq, olinucs, length, k, prop, supInfo)
    while j < limit or j == limit:
        current = getSpecificValue(seq[j - 1], olinucs, prop, supInfo)
        # hn(seq[j-1],prop,SupFileName)
        partOne = current - avg
        next = getSpecificValue(seq[j + l - 1], olinucs, prop, supInfo)
        # hn(seq[j+l-1],prop,SupFileName)
        partTwo = next - avg
        top = top + (partOne * partTwo)
        j = j + 1
    top = top / limit
    limit2 = length - k + 1
    bottom = 0
    b = 1
    while b < limit2 or b == limit2:
        current = getSpecificValue(seq[b - 1], olinucs, prop, supInfo)
        # hn(seq[b-1],prop,SupFileName)
        bottom = bottom + ((current - avg) * (current - avg))
        b = b + 1
    bottom = bottom / limit2
    final = top / bottom
    return final

def autocorrelation(autoc, inputfile, props, k, l, alphabet):
    if not props:
        error_info = 'Error, The phyche_list, extra_index_file and all_prop can\'t be all False.'
        raise ValueError(error_info)

    input_data = open(inputfile, 'r')
    sequences =get_data(input_data,alphabet)
    # Getting supporting info from files
    if k == 2 and alphabet == index_list.RNA:
        SupFileName = './data/Supporting_Information_S1_RNA.txt'
    elif k == 2 and alphabet == index_list.DNA:
        SupFileName = './data/Supporting_Information_S1_DNA.txt'
    elif k == 3 and alphabet == index_list.DNA:
        SupFileName = './data/Supporting_Information_S3_DNA.txt'
    SupFile = open(SupFileName, 'r')
    supInfo = SupFile.read()
    o = re.search('Physicochemical properties\,(.+)\n', supInfo)
    olinucs = ''
    if o:
        olinucs = o.group(1).rstrip()
    SupFile.close()
    # Writing to output file
    m = 0
    vectors = []
    for sequence in sequences:
        length = len(sequence)
        seq = sepSequence(sequence, k)
        values = []
        for prop in props:
            if autoc.upper() == 'MAC':
                value = float("%.3f" % moran(seq, olinucs, length, k, l, prop, supInfo))
                values.append(value)
            elif autoc.upper() == 'GAC':
                value = float("%.3f" % geary(seq, olinucs, length, k, l, prop, supInfo))
                values.append(value)
            elif autoc.upper() == 'NMBAC':
                value = float("%.3f" % moreau(seq, olinucs, length, k, l, prop, supInfo))
                values.append(value)
        vectors.append(values)
        m += 1
    return vectors


#====================================================================================================

def acc(input_data, k, lag, phyche_list, alphabet, extra_index_file=None, all_prop=False, theta_type=1):
    """This is a complete acc in PseKNC.

    :param k: int, the value of k-tuple.
    :param phyche_list: list, the input physicochemical properties list.
    :param extra_index_file: a file path includes the user-defined phyche_index.
    :param all_prop: bool, choose all physicochemical properties or not.
    :param theta_type: the value 1, 2 and 3 for ac, cc or acc.
    """
    phyche_list = get_phyche_list(k, phyche_list,
                                  extra_index_file=extra_index_file, alphabet=alphabet, all_prop=all_prop)
    # print(phyche_list)
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
        # print(phyche_vals)
        if extra_index_file is not None:
            phyche_vals.extend(extend_aaindex(extra_index_file))

    seqs = get_data(input_data, alphabet)
    if alphabet == index_list.PROTEIN:
        # Transform the data format to dict {acid: [phyche_vals]}.
        phyche_keys = phyche_vals[0].index_dict.keys()
        phyche_vals = [e.index_dict.values() for e in phyche_vals]
        new_phyche_vals = zip(*[e for e in phyche_vals])
        phyche_vals = {key: list(val) for key, val in zip(phyche_keys, new_phyche_vals)}

    if theta_type == 1:
        return make_ac_vec(seqs, lag, phyche_vals, k)
    elif theta_type == 2:
        return make_cc_vec(seqs, lag, phyche_vals, k)
    elif theta_type == 3:
        return make_acc_vec(seqs, lag, phyche_vals, k)


def make_ac_vec(sequence_list, lag, phyche_value, k):

    # Get the length of phyche_vals.
    phyche_values = list(phyche_value.values())
    len_phyche_value = len(phyche_values[0])

    vec_ac = []
    for sequence in sequence_list:
        len_seq = len(sequence)
        each_vec = []

        for temp_lag in range(1, lag + 1):
            for j in range(len_phyche_value):

                # Calculate average phyche_value for a nucleotide.
                ave_phyche_value = 0.0
                for i in range(len_seq - temp_lag - k + 1):
                    nucleotide = sequence[i: i + k]
                    ave_phyche_value += float(phyche_value[nucleotide][j])
                ave_phyche_value /= len_seq

                # Calculate the vector.
                temp_sum = 0.0
                for i in range(len_seq - temp_lag - k + 1):
                    nucleotide1 = sequence[i: i + k]
                    nucleotide2 = sequence[i + temp_lag: i + temp_lag + k]
                    temp_sum += (float(phyche_value[nucleotide1][j]) - ave_phyche_value) * (
                        float(phyche_value[nucleotide2][j]))

                each_vec.append(round(temp_sum / (len_seq - temp_lag - k + 1), 8))
        vec_ac.append(each_vec)

    return vec_ac


def make_cc_vec(sequence_list, lag, phyche_value, k):
    phyche_values = list(phyche_value.values())
    len_phyche_value = len(phyche_values[0])

    vec_cc = []
    for sequence in sequence_list:
        len_seq = len(sequence)
        each_vec = []

        for temp_lag in range(1, lag + 1):
            for i1 in range(len_phyche_value):
                for i2 in range(len_phyche_value):
                    if i1 != i2:
                        # Calculate average phyche_value for a nucleotide.
                        ave_phyche_value1 = 0.0
                        ave_phyche_value2 = 0.0
                        for j in range(len_seq - temp_lag - k + 1):
                            nucleotide = sequence[j: j + k]
                            ave_phyche_value1 += float(phyche_value[nucleotide][i1])
                            ave_phyche_value2 += float(phyche_value[nucleotide][i2])
                        ave_phyche_value1 /= len_seq
                        ave_phyche_value2 /= len_seq

                        # Calculate the vector.
                        temp_sum = 0.0
                        for j in range(len_seq - temp_lag - k + 1):
                            nucleotide1 = sequence[j: j + k]
                            nucleotide2 = sequence[j + temp_lag: j + temp_lag + k]
                            temp_sum += (float(phyche_value[nucleotide1][i1]) - ave_phyche_value1) * \
                                        (float(phyche_value[nucleotide2][i2]) - ave_phyche_value2)
                        each_vec.append(round(temp_sum / (len_seq - temp_lag - k + 1), 8))

        vec_cc.append(each_vec)

    return vec_cc


def make_acc_vec(seqs, lag, phyche_values, k):
    from functools import reduce
    zipped = list(zip(make_ac_vec(seqs, lag, phyche_values, k), make_cc_vec(seqs, lag, phyche_values, k)))
    return [reduce(lambda x, y: x + y, e) for e in zipped]


#--------------------------------------------------------------------------
#PDT method
#--------------------------------------------------------------------------

def pdt_cmd(inputfile, lamada):
    """Concatenation of pdt command.
    :param inputfile: the input sequence file in FASTA format.
    :param lamada: the value of parameter lamada.
    """
    if sys.platform.startswith('win'):
        pdt_cmd = '.\\pdt\\pdt.exe'
    else:
        pdt_cmd = './pdt/pdt'
        #os.chmod(pdt_cmd, stat.S_IRWXU)
        os.chmod(pdt_cmd, 0777)

    input_file = inputfile
    aaindex_file = './pdt/aaindex_norm.txt'
    lamada_param = str(lamada)
    output_file = ''.join([os.path.splitext(input_file)[0], '_temp', os.path.splitext(input_file)[1]])

    cmd = ''.join([pdt_cmd, ' ', input_file, ' ', aaindex_file, ' ', lamada_param, ' ', output_file])
    #print cmd
    subprocess.call(cmd, shell=True)
    return os.path.abspath(output_file)


def pdt(inputfile, lamada):
    """Execute pdt command and generate feature vectors.
    :param inputfile: the input sequence file in FASTA format.
    :param lamada: the value of parameter lamada.
    """
    #inputfile = os.path.abspath(inputfile)
    #cwd = os.getcwd()
    #os.chdir('./PDT')
    output_file = pdt_cmd(inputfile, lamada)
    #os.chdir(cwd)

    vector_list = []
    with open(output_file, 'r') as f:
        for line in f:
            temp_list = line.strip().split('\t')
            vector = [round(float(elem), 3) for elem in temp_list]
            vector_list.append(vector)
    #print vector_list
    return vector_list




def main(args):
    """The main process of autocorrelation methods.
    :param args: an object of the arguments.
    """
    file_list = args.inputfiles

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


    if args.method.upper() not in ['MAC', 'GAC', 'NMBAC', 'PDT']:
        for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
            with open(input_file) as f:
                k = read_k(args.alphabet, args.method, 0)

            # Get index_list.
                if args.i is not None:
                    from pse import read_index
                    ind_list = read_index(args.i)
                else:
                    ind_list = []

                default_e = []
            # Set Pse default index_list.
                if args.alphabet == 'DNA':
                    alphabet_list = index_list.DNA
                    if k == 2:
                        default_e = const.DI_INDS_6_DNA
                    elif k == 3:
                        default_e = const.TRI_INDS_DNA
                elif args.alphabet == 'RNA':
                    alphabet_list = index_list.RNA
                    default_e = const.DI_INDS_RNA
                elif args.alphabet == 'Protein':
                    alphabet_list = index_list.PROTEIN
                    default_e = const.INDS_3_PROTEIN
                else:
                    print 'The alphabet should be DNA, RNA or Protein.'
                    return False

                theta_type = 1
                if args.method in const.METHODS_AC:
                    theta_type = 1
                elif args.method in const.METHODS_CC:
                    theta_type = 2
                elif args.method in const.METHODS_ACC:
                    theta_type = 3
                else:
                    print("Method error!")


            # ACC.

                if args.e is None and len(ind_list) == 0 and args.a is False:
                # Default Pse.
                    res = acc(f, k, args.lag, default_e, alphabet_list,
                              extra_index_file=args.e, all_prop=args.a, theta_type=theta_type)
                else:
                    res = acc(f, k, args.lag, ind_list, alphabet_list,
                              extra_index_file=args.e, all_prop=args.a, theta_type=theta_type)
            write_to_file(res, output_format, label, output_file)

    if args.method.upper() in ['MAC', 'GAC', 'NMBAC']:
        if args.lamada < 0 or args.lamada > 10:
            print 'The value of lamada should be larger than 0 and smaller than 10.'
            return False
        if args.a is None:
            args.a == False
        elif args.alphabet == 'DNA':
            args.alphabet = index_list.DNA
            if args.oli == 0:
                if args.a == True:
                    for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
                        res = autocorrelation(autoc=args.method, inputfile=input_file, props=const.ALL_DI_DNA_IND, k=2, l=args.lamada, alphabet=args.alphabet)
                        write_to_file(res, output_format, label, output_file)
                elif args.a == False:
                    for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
                        res = autocorrelation(autoc=args.method, inputfile=input_file, props=const.DEFAULT_DI_DNA_IND, k=2, l=args.lamada, alphabet=args.alphabet)
                        write_to_file(res, output_format, label, output_file)
            if args.oli == 1:
                if args.a == True:
                    for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
                        res = autocorrelation(autoc=args.method, inputfile=input_file, props=const.ALL_TRI_DNA_IND, k=3, l=args.lamada, alphabet=args.alphabet)
                        write_to_file(res, output_format, label, output_file)
                elif args.a == False:
                    for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
                        res = autocorrelation(autoc=args.method, inputfile=input_file, props=const.DEFAULT_TRI_DNA_IND, k=3, l=args.lamada, alphabet=args.alphabet)
                        write_to_file(res, output_format, label, output_file)
        elif args.alphabet == 'RNA':
            args.alphabet = index_list.RNA
            if args.a == True:
                for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
                    res = autocorrelation(autoc=args.method, inputfile=input_file, props=const.ALL_RNA_IND, k=2, l=args.lamada, alphabet=args.alphabet)
                    write_to_file(res, output_format, label, output_file)
            elif args.a == False:
                for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
                    res = autocorrelation(autoc=args.method, inputfile=input_file, props=const.DEFAULT_RNA_IND, k=2, l=args.lamada, alphabet=args.alphabet)
                    write_to_file(res, output_format, label, output_file)

    if args.method.upper() == 'PDT':
        if args.alphabet != 'Protein':
            print 'PDT method is only available for Protein sequences.'
            return False
        else:
            if args.lamada < 1 or args.lamada > 15:
                print 'The value of -lamada should be larger than 0 and smaller than 16.'
                return False
            else:
                for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
                    res = pdt(input_file, args.lamada)
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

    parse = argparse.ArgumentParser(description="This is acc module for generate acc vector.",
                                    formatter_class=RawTextHelpFormatter)

    parse.add_argument('inputfiles', nargs='*',
                       help="The input files in FASTA format.")
    parse.add_argument('-out', nargs='*',
                       help="The output files for storing feature vectors.")
    parse.add_argument('alphabet', choices=['DNA', 'RNA', 'Protein'],
                       help="The alphabet of sequences.")
    parse.add_argument('method', type=str,
                       help="The method name of autocorrelation.")

    parse.add_argument('-lag', type=int, default=2,
                       help="The value of lag.")
    parse.add_argument('-lamada', type=int, default=1,
                       help="The value of lamada. default=1")
    parse.add_argument('-oli', type=int, default=0, choices=[0,1],
                       help="Choose one kind of Oligonucleotide: 0 represents dinucleotide;\n"
                       "1 represents trinucleotide.")
    parse.add_argument('-i',
                       help="The indices file user choose.\n"
                            "Default indices:\n"
                            "DNA dinucleotide: Rise, Roll, Shift, Slide, Tilt, Twist.\n"
                            "DNA trinucleotide: Dnase I, Bendability (DNAse).\n"
                            "RNA: Rise, Roll, Shift, Slide, Tilt, Twist.\n"
                            "Protein: Hydrophobicity, Hydrophilicity, Mass.")
    parse.add_argument('-e',
                       help="The user-defined indices file.")
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
                       "For multiclass classification problem, the labels can be set as a list integers.")


    args = parse.parse_args()
    # print(args)

    if check_args(args, 'acc.py'):
        print("Calculating...")
        start_time = time.time()
        main(args)
        print("Done.")
        print("Used time: %.2fs" % (time.time() - start_time))

