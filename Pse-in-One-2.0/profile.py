#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys, os, subprocess
import heapq
import time
import cPickle
import random
import multiprocessing as mul
from util import Seq, write_to_file
from nac import make_kmer_list, dr_method
from acc import pdt
import const


def is_fasta_and_protein(seq):
    """Judge if the seq is in fasta format and protein sequences.
    :param seq: Seq object
    Return True or False.
    """
    if not seq.name:
        error_info = 'Error, sequence ' + str(seq.no) + ' has no sequence name.'
        print(seq)
        sys.stderr.write(error_info)
        return False
    # if -1 != seq.name.find('>'):
    #     error_info = 'Error, sequence ' + str(seq.no) + ' name has > character.'
    #     sys.stderr.write(error_info)
    #    return False
    if 0 == seq.length:
        error_info = 'Error, sequence ' + str(seq.no) + ' is null.'
        sys.stderr.write(error_info)
        return False
    for elem in seq.seq:
        if elem not in const.PROTEIN and elem !='x':
            error_info = 'Sorry, sequence ' + str(seq.no) \
                         + ' has character ' + str(elem) + '.(The character must be ' + const.PROTEIN + ').'
            sys.stderr.write(error_info)
            return False
    return True


def check_and_save(filename):
    """Read the input file and store as Seq objects.
    :param filename: the input protein sequence file.
    return  an iterator.
    """
    name, seq = '', ''
    count = 0
    seq_list = []
    with open(filename) as f:
    #lines = f.readlines()
        for line in f:
            if not line:
                break

            if '>' == line[0]:
                if 0 != count or (0 == count and seq != ''):
                    if is_fasta_and_protein(Seq(name, seq, count)):
                        yield Seq(name, seq, count)
                    else:
                        sys.exit(0)

                seq = ''
                name = line[1:].strip()
                count += 1
            else:
                if 'x' in line or 'X' in line:
                    line = line.replace('x', '')
                    line = line.replace('X', '')
                seq += line.strip()

        #count += 1
        if is_fasta_and_protein(Seq(name, seq, count)):
            yield Seq(name, seq, count)
        else:
            sys.exit(0)




def sep_file(filename):
    """seperate the input file. One sequence in one file.
    :param filename: the input file.
    """
    dirname, suff  = os.path.splitext(filename)
    #dirname = ''.join([const.TEMP_DIR, os.path.split(dirname)[1]])
    try:
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        else:
            randstr = str(random.randint(0, 99999))
            dirname = dirname + '_' + randstr
            os.mkdir(dirname)
    except OSError:
        dir_path = os.path.split(dirname)[0]
        new_dir = str(random.randint(10000, 99999))
        dirname = ''.join([dir_path, '/',  new_dir])
        if not os.path.isdir(dirname):
            os.mkdir(dirname)

    seq_name = []
    for seq in check_and_save(filename):
        #print seq.no
        seq_name.append(seq.name)
        seq_file = dirname + '/' + str(seq.no) + '.txt'
        with open(seq_file, 'w') as f:
            f.write('>')
            f.write(str(seq.name))
            f.write('\n')
            f.write(str(seq.seq))
    return os.path.abspath(dirname), seq_name



def produce_one_frequency(fasta_file, xml_file, pssm_file, semph):
    """Produce fequency profile for one sequence using psiblast.
    :param fasta_file: the file storing one sequence.
    :param xml_file: the generated xml file by psiblast.
    :param pssm_file: the generated pssm file by psiblast.
    :param semph: the semaphore used for multiprocessing.
    """
    #print type(os.listdir(dirname))
    semph.acquire()
    if sys.platform.startswith('win'):
        psiblast_cmd = '.\\psiblast\\psiblast.exe'
        #psiblast_cmd = 'psiblast.exe'
    else:
        psiblast_cmd = './psiblast/psiblast'
        os.chmod(psiblast_cmd, 0777)
    BLAST_DB = './psiblast/nrdb90/nrdb90'
    #BLAST_DB = 'nrdb90/nrdb90'
    #BLAST_DB = 'nr/nr'
    #BLAST_DB = 'uniref50/uniref50'

    #xml_file = 'test.xml'
    #evalue_threshold = 1
    evalue_threshold = 0.001
    num_iter = 3
    #pssm_file = 'test.pssm'
    outfmt_type = 5

    cmd = ' '.join([psiblast_cmd,
                   '-query ' + fasta_file,
                   '-db ' + BLAST_DB,
                    '-out ' + xml_file,
                   '-evalue ' + str(evalue_threshold),
                    '-num_iterations ' + str(num_iter),
                   '-num_threads ' + '5',
                    '-out_ascii_pssm ' + pssm_file,
                    '-outfmt ' + str(outfmt_type)
                    ]
                )
    #print cmd

    #cwd = os.getcwd()
    #os.chdir('./psiblast')
    subprocess.call(cmd, shell=True)
    semph.release()
    #os.chdir(cwd)


def produce_all_frequency(dirname, process_num):
    """Produce frequency profile for all the sequences.
    :param dirname: the directory used to store the generated files.
    :param process_num: the number of processes used for multiprocessing.
    """
    sequence_files = []
    for i in os.listdir(dirname):
        seq_full_path = ''.join([dirname, '/', i])
        if os.path.isfile(seq_full_path):
            sequence_files.append(seq_full_path)
    #print sequence_files

    process_list = []
    semph = mul.Semaphore(process_num)

    xml_dir = ''.join([dirname, '/xml'])
    if not os.path.isdir(xml_dir):
        os.mkdir(xml_dir)
    pssm_dir = ''.join([dirname, '/pssm'])
    if not os.path.isdir(pssm_dir):
        os.mkdir(pssm_dir)
    for seq_file in sequence_files:
        #fasta_file = seq_file
        name = os.path.splitext(os.path.split(seq_file)[1])[0]
        xml_file = ''.join([xml_dir, '/', name, '.xml'])
        pssm_file = ''.join([pssm_dir, '/', name, '.pssm'])
        process_list.append(mul.Process(target=produce_one_frequency,
                                        args=(seq_file, xml_file, pssm_file, semph)))
    #cwd = os.getcwd()
    #os.chdir('./psiblast')
    for process in process_list:
        process.start()
    for process in process_list:
        process.join()
    #os.chdir(cwd)
    return pssm_dir
        #print xml_file
        #print pssm_file

def produce_one_top_n_gram(pssm_file, n):
    """Produce top-n-gram for one pssm file.
    :param pssm_file: the pssm file used to generate top-n-gram.
    :param n: the n most frequent amino acids in the amino acid frequency profiles.
    """
    tng_list = []
    new_alpha_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
        'P', 'S', 'T', 'W', 'Y', 'V']
    with open(pssm_file, 'r') as f:
        count = 0
        for line in f:
            count += 1
            if count <= 3 or len(line.strip().split()) < 20:
                pass
            else:
                line = line.strip().split()[22:-2]
                line = map(eval, line)
                #print line
                data = heapq.nlargest(n, enumerate(line), key=lambda x:x[1])
                #print data
                if n == 1:
                    new_alpha = new_alpha_list[data[0][0]]
                else:
                    new_alpha = ''
                    indices, val = zip(*data)
                    for i in indices:
                        new_alpha += new_alpha_list[i]
                tng_list.append(new_alpha)

    #print len(tng_list)
    return tng_list


def produce_tng_blosum(seq_file, blosum_dict, n):
    """Generate top-n-gram by blosum62 matrix.
    :param seq_file: the sequence file containing one sequence.
    :param blosum_dict: the dict which stores the blosum62 matrix.
    :param n: the n most frequent amino acids in the amino acid frequency profiles.
    """
    tng_list = []
    with open(seq_file, 'r') as f:
        for line in f:
            if line.strip().startswith('>'):
                continue
            else:
                line = line.strip()
                for amino in line:
                    amino_blosum = blosum_dict[amino]
                    data = heapq.nlargest(n, enumerate(amino_blosum), key=lambda x:x[1])
                    if n == 1:
                        index = data[0][0]
                        new_alpha = blosum_dict['alphas'][index]
                    else:
                        new_alpha = ''
                        indices, val = zip(*data)
                        for i in indices:
                            new_alpha += blosum_dict['alphas'][i]
                    tng_list.append(new_alpha)
    return tng_list


def produce_top_n_gram(pssm_dir, seq_name, n):
    """Produce top-n-gram for all the pssm files.
    :param pssm_dir: the directory used to store pssm files.
    :param seq_name: the name of sequences.
    :param n: the n most frequent amino acids in the amino acid frequency profiles.
    """
    dirname = os.path.split(pssm_dir)[0]
    fasta_name = os.path.split(dirname)[1]
    final_result = ''.join([dirname, '/final_result'])
    if not os.path.isdir(final_result):
            os.mkdir(final_result)
    pssm_files = []
    dir_list = os.listdir(pssm_dir)
    index_list = []
    for elem in dir_list:
        pssm_full_path = ''.join([pssm_dir, '/', elem])
        name, suffix = os.path.splitext(elem)
        if os.path.isfile(pssm_full_path) and suffix == '.pssm':
            index_list.append(int(name))

    index_list.sort()

    #return
    if len(index_list) != len(seq_name):
        with open('./psiblast/blosum62.pkl', 'rb') as f:
            blosum_dict = cPickle.load(f)
    tng_all_list = []

    #print tng_all_list
    for i in xrange(1, len(seq_name)+1):
        if i in index_list:
            pssm_full_path = ''.join([pssm_dir, '/', str(i), '.pssm'])
            tng = produce_one_top_n_gram(pssm_full_path, n)
        else:
            seq_file = ''.join([dirname, '/', str(i), '.txt'])
            tng = produce_tng_blosum(seq_file, blosum_dict, n)
        tng_all_list.append(tng)

    tng_file_name = ''.join([final_result, '/', fasta_name, '_new.txt'])
    with open(tng_file_name, 'w') as f:
        for index, tng in enumerate(tng_all_list):
            f.write('>')
            f.write(seq_name[index])
            f.write('\n')
            for elem in tng:
                f.write(elem)
                f.write(' ')
            f.write('\n')
    #print tng_file_name
    return tng_all_list




def top_n_gram(inputfile, n, process_num):
    """Generate top-n-gram list.
    :param inputfile: input sequence file in FASTA format.
    :param n: the n most frequent amino acids in the amino acid frequency profiles.
    :param process_num: the number of processes used for multiprocessing.
    """
    dirname, seq_name = sep_file(inputfile)
    pssm_dir = produce_all_frequency(dirname, process_num)
    tng_all_list = produce_top_n_gram(pssm_dir, seq_name, n)

    return tng_all_list, seq_name


def gen_vector(tng_list, n):
    """Generate feature vectors based on the top-n-gram.
    :param tng_list: the generated top-n-gram list.
    :param n: the n most frequent amino acids in the amino acid frequency profiles.
    """
    gram_list = make_kmer_list(n, const.PROTEIN)
    #print tng_list
    #print gram_list
    vector_list = []
    for tng in tng_list:
        vec_len = len(tng)
        #print vec_len
        vector = []
        for elem in gram_list:
            gram_count = tng.count(elem)
            occur_freq = round((gram_count * 1.0) / vec_len, 4)
            vector.append(occur_freq)
        vector_list.append(vector)

    return vector_list

#-------------------------------------------------------------------------------------
# PDT-Profile start
#-------------------------------------------------------------------------------------

def convert_tng_to_fasta(tng_list, seq_name, origin_file_name):
    """Convert top-n-gram to fasta format.
    :param tng_list: the generated top-n-gram list.
    :param seq_name: the name of sequences.
    :param origin_file_name: the name of the input file in FASTA format.
    """
    name = os.path.split(origin_file_name)[1]
    name = os.path.splitext(name)
    tng_fasta = ''.join([const.TEMP_DIR, name[0], '_tng.txt'])
    with open(tng_fasta, 'w') as f:
        for index, tng in enumerate(tng_list):
            f.write('>')
            f.write(seq_name[index])
            f.write('\n')
            for elem in tng:
                f.write(elem)
            f.write('\n')
    return tng_fasta


def pdt_profile(inputfile, n, lamada, process_num):
    """Generate PDT-Profile features.
    :param inputfile: input sequence file in FASTA format.
    :param n: the n most frequent amino acids in the amino acid frequency profiles.
    :param lamada: the distance between two amino acids.
    :param process_num: the number of processes used for multiprocessing.
    """
    tng_list, seq_name = top_n_gram(inputfile, n, process_num)
    tng_fasta = convert_tng_to_fasta(tng_list, seq_name, inputfile)
    return pdt(tng_fasta, lamada)

#-------------------------------------------------------------------------------------
# PDT-Profile end
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# DT start
#-------------------------------------------------------------------------------------

def dt_method(inputfile, max_dis, process_num):
    """Generate DT method feature vectors.
    :param inputfile: input sequence file in FASTA format.
    :param max_dis: the maximum distance between top-1-gram pairs.
    :param process_num: the number of processes used for multiprocessing.
    """
    tng_list, seq_name = top_n_gram(inputfile, 1, process_num)
    tng_fasta = convert_tng_to_fasta(tng_list, seq_name, inputfile)
    return dr_method(tng_fasta, max_dis)

#-------------------------------------------------------------------------------------
# DT end
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# ACC-PSSM, AC-PSSM, CC-PSSM start
#-------------------------------------------------------------------------------------

def blosum_pssm(seq_file, new_blosum_dict, dirname):
    """Generate pssm file using blosum62 matrix.
    :param seq_file: the sequence file containing one sequence.
    :param new_blosum_dict: the blosum62 dict after processing.
    :param dirname: the directory name for storing the generated files.
    """
    pssm_list = []
    with open(seq_file, 'r') as f:
        for line in f:
            if line.strip().startswith('>'):
                continue
            else:
                for index, i in enumerate(line.strip()):
                    blosum_list = map(str, new_blosum_dict[i])
                    blosum_str = ' '.join(blosum_list)
                    line_str = ' '.join([str(index+1), i, blosum_str])
                    pssm_list.append(line_str)

    blosum_dir = ''.join([dirname, '/blosum_pssm'])
    if not os.path.isdir(blosum_dir):
        os.mkdir(blosum_dir)
    seq_file_name = os.path.splitext(seq_file)[0]
    seq_file_name = os.path.split(seq_file_name)[1]
    blosum_file = ''.join([blosum_dir, '/', seq_file_name, '.pssm'])
    with open(blosum_file, 'w') as f:
        for i in range(3):
            f.write('\n')
        for line in pssm_list:
            f.write(line)
            f.write('\n')
    return os.path.abspath(blosum_file)



def read_blosum():
    """Read blosum dict and delete some keys and values."""
    with open('./psiblast/blosum62.pkl', 'rb') as f:
        blosum_dict = cPickle.load(f)

    temp = blosum_dict.pop('*')
    temp = blosum_dict.pop('B')
    temp = blosum_dict.pop('Z')
    temp = blosum_dict.pop('X')
    temp = blosum_dict.pop('alphas')

    for key in blosum_dict:
        for i in range(4):
            temp = blosum_dict[key].pop()
    return blosum_dict



def acc_pssm_cmd(pssm_file, lag, acc_out_file):
    """ACC-PSSM command.
    :param pssm_file: the .pssm file.
    :param lag: the distance between two amino acids.
    :param acc_out_file: the output file of the acc program.
    """
    if sys.platform.startswith('win'):
        acc_cmd = '.\\acc_pssm\\acc.exe'
    else:
        acc_cmd = './acc_pssm/acc'
        os.chmod(acc_cmd, 0777)

    cmd = ' '.join([acc_cmd, ' ', str(lag), ' ', pssm_file, ' ', acc_out_file])
    subprocess.call(cmd, shell=True)



def sep_acc_vector(acc_out_file):
    """Seperate acc_out_file and output the acc, ac, cc vectors.
    :param acc_out_file: the output file of the acc program.
    """
    acc_vec_list = []
    ac_vec_list = []
    cc_vec_list = []
    with open(acc_out_file, 'r') as f:
        for line in f:
            line = round(float(line.strip()), 3)
            #line = float(line.strip())
            acc_vec_list.append(line)

    for i in range(0, len(acc_vec_list), 400):
        ac_vec_list.extend(acc_vec_list[i:i+20])
        cc_vec_list.extend(acc_vec_list[i+20:i+400])

    return acc_vec_list, ac_vec_list, cc_vec_list



def make_acc_pssm_vector(inputfile, lag, vec_type, process_num):
    """Generate ACC, AC, CC feature vectors.
    :param inputfile: input sequence file in FASTA format.
    :param lag: the distance between two amino acids.
    :param vec_type: the type of the vectors generated, ACC-PSSM, AC-PSSM
    or CC-PSSM.
    :param process_num: the number of processes used for multiprocessing.
    """
    dirname, seq_name = sep_file(inputfile)
    pssm_dir = produce_all_frequency(dirname, process_num)

    dir_list = os.listdir(pssm_dir)
    index_list = []
    for elem in dir_list:
        pssm_full_path = ''.join([pssm_dir, '/', elem])
        name, suffix = os.path.splitext(elem)
        if os.path.isfile(pssm_full_path) and suffix == '.pssm':
            index_list.append(int(name))

    index_list.sort()
    if len(index_list) != len(seq_name):
        new_blosum_dict = read_blosum()

    acc_out_fold = dirname + '/acc_out'
    acc_vectors = []
    ac_vectors = []
    cc_vectors = []
    if not os.path.isdir(acc_out_fold):
        os.mkdir(acc_out_fold)
    for i in xrange(1, len(seq_name)+1):
        if i in index_list:
            pssm_full_path = ''.join([pssm_dir, '/', str(i), '.pssm'])
        else:
            seq_file = ''.join([dirname, '/', str(i), '.txt'])
            pssm_full_path = blosum_pssm(seq_file, new_blosum_dict, dirname)
        acc_out_file = ''.join([acc_out_fold, '/', str(i), '.out'])
        acc_pssm_cmd(pssm_full_path, lag, acc_out_file)
        acc_vec_list, ac_vec_list, cc_vec_list = sep_acc_vector(acc_out_file)
        acc_vectors.append(acc_vec_list)
        ac_vectors.append(ac_vec_list)
        cc_vectors.append(cc_vec_list)
    if vec_type == 'acc':
        return acc_vectors
    elif vec_type == 'ac':
        return ac_vectors
    elif vec_type == 'cc':
        return cc_vectors
    else:
        return False

#-------------------------------------------------------------------------------------
# ACC-PSSM, AC-PSSM, CC-PSSM end
#-------------------------------------------------------------------------------------
def main(args):
    """The main process of profile-based features.
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

    cpu_core = mul.cpu_count()
    if args.cpu is None:
        process_num = 1
    elif 0 < args.cpu <= cpu_core:
        process_num = args.cpu
    elif args.cpu < 0 or args.cpu > cpu_core:
        print 'Error: The value of -cpu should be larger than 0'
        print 'and less than or equal to the number of cpu core in your computer.'
        return False

    if args.method.upper() == 'TOP-N-GRAM':
        for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
            tng_list = top_n_gram(input_file, args.n, process_num)[0]
            res = gen_vector(tng_list, args.n)
            write_to_file(res, output_format, label, output_file)
    elif args.method.upper() == 'PDT-PROFILE':
        if args.lamada < 1:
            print 'The value of -lamada should be larger than 0.'
            return False
        else:
            for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
                res = pdt_profile(input_file, args.n, args.lamada, process_num)
                write_to_file(res, output_format, label, output_file)
    elif args.method.upper() == 'DT':
        if args.max_dis < 1:
            print 'The value of -max_dis should be larger than 0.'
            return False
        else:
            for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
                res = dt_method(input_file, args.max_dis, process_num)
                write_to_file(res, output_format, label, output_file)
    elif args.method.upper() in ['ACC-PSSM', 'AC-PSSM', 'CC-PSSM']:
        if args.method.upper() == 'ACC-PSSM':
            vec_type = 'acc'
        elif args.method.upper() == 'AC-PSSM':
            vec_type = 'ac'
        elif args.method.upper() == 'CC-PSSM':
            vec_type = 'cc'
        if args.lag < 1:
            print 'The value of -lag should be larger than 0.'
            return False
        else:
            for input_file, output_file, label in zip(file_list, outputfile_list, label_list):
                res = make_acc_pssm_vector(input_file, args.lag, vec_type, process_num)
                write_to_file(res, output_format, label, output_file)

    else:
        print 'Method error!'
        return False
    if len(outputfile_list) != 0:
        for index, output_file in enumerate(outputfile_list):
            out_with_full_path = os.path.abspath(output_file)
            if os.path.isfile(out_with_full_path):
                if index == 0:
                    print 'The output file(s) can be found here:'
                print out_with_full_path






if __name__ == '__main__':
    # for i in check_and_save('protein_412.txt'):
    #    print i.name
    #    print i.seq
    #    print i.no
    # dirname = sep_file('protein_test.txt')
    # #print dirname
    # res = produce_all_frequency(dirname)
    # print res
    #produce_frequency('1.txt')
    #produce_one_top_n_gram('1.pssm', 2)

    #vec = top_n_gram('protein_test.txt', 2)
    #write_to_file(vec, 'svm', '+1', 'vectors.txt')

    import argparse
    from argparse import RawTextHelpFormatter

    parse = argparse.ArgumentParser(description="This is profile module for generate profile-based vectors.",
                                    formatter_class=RawTextHelpFormatter)
    parse.add_argument('inputfiles', nargs='*',
                       help="The input files in FASTA format.")
    parse.add_argument('-out', nargs='*',
                       help="The output files for storing feature vectors.")
    parse.add_argument('method', type=str,
                        help="The method name of profile-based features. {Top-n-gram,PDT-Profile,DT,AC-PSSM,CC-PSSM,ACC-PSSM}")
    parse.add_argument('-n', type=int, choices=[1, 2, 3],
                       help="For Top-n-gram,PDT-Profile methods. The value of top-n-gram.")
    parse.add_argument('-lamada', type=int, default=1,
                       help="The value of lamada. default=1")
    parse.add_argument('-max_dis', type=int, default=3,
                       help="For DT methods. The max distance value of residues.")
    parse.add_argument('-lag', type=int, default=2,
                       help="For ACC-PSSM, AC-PSSM and CC-PSSM methods. The value of lag.")
    parse.add_argument('-f', default='tab', choices=['tab', 'svm', 'csv'],
                       help="The output format (default = tab).\n"
                            "tab -- Simple format, delimited by TAB.\n"
                            "svm -- The libSVM training data format.\n"
                            "csv -- The format that can be loaded into a spreadsheet program.")
    parse.add_argument('-labels', nargs='*',
                       help="The labels of the input files.\n"
                       "For binary classification problem, the labels can only be '+1' or '-1'.\n"
                       "For multiclass classification problem, the labels can be set as a list of integers.")
    parse.add_argument('-cpu', type=int,
                       help="The maximum number of CPU core used for multiprocessing in\n"
                       "generating frequency profile. Default is 1.")

    args = parse.parse_args()
    print("Calculating...")
    start_time = time.time()
    main(args)
    print("Done.")
    print("Used time: %.2fs" % (time.time() - start_time))