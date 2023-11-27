# CircDC
1. Description
 Circular RNAs (circRNAs) can regulate microRNA activity and are related to various diseases, such as cancer. Functional research on circRNAs is the focus of scientific research. Accurate identification of circRNAs is important for gaining insight into their functions.


We developed a novel framework, CircDC, for classifying circRNAs from other lncRNAs. CircDC uses four different feature encoding schemes and adopts a multilayer convolutional neural network and bidirectional long short-term memory network to learn high-order feature representation and make circRNA predictions. 

2. Availability

2.1 Datasets and source code are available at:
https://github.com/nmt315320/CircDC.git.
Data：circRNA data, .bed format,circRNA data, bed format
      hg38.fa---Because the data is relatively large, you need to download it yourself.

2.1 Local running

2.1.1 Environment
Before running, please make sure the following packages are installed in Python environment:

gensim==3.4.0
pysam

pigwig

pandas==1.1.3

tensorflow==2.3.0

python==3.7.3

biopython==1.7.8

numpy==1.19.2

For convenience, we strongly recommended users to install the Anaconda Python 3.7.3 (or above) in your local computer.

2.1.2 Additional requirements
One additional file, namely hg.38.fa, is needed for CircDC, we did not provide this in the source code packages because of the license restriction. This file can be acquired at the following links:

hg38: wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz

2.1.3 Running

Changing working dir to CircDC-master, and then running the following command:

python CircDC.py -i testing.fasta -o results.csv 

-i: name of input_file in fasta format # folder “sequence” is the default file path of the input_file

-o name of output_file # folder “results” is the default file path for result save.

The optional parameters are:
    parser.add_argument('--data_dir', type=str, default='data/', metavar='<data_directory>',
                        help='Under this directory, you will have descriptors files ')

    parser.add_argument('--train', type=bool, default=True, help='use this option for training model')

    parser.add_argument('--model_dir', type=str, default='models/',
                        help='The directory to save the trained models for future prediction')

    parser.add_argument('--predict', type=bool, default=False,
                        help='Predicting circular RNAs. if using train, then it will be False')

    parser.add_argument('--out_file', type=str, default='prediction.txt',
                        help='The output file used to store the prediction probability of testing data')

    parser.add_argument('--rcm', type=bool, default=True, help='The modularity of RCM')

    parser.add_argument('--cons', type=bool, default=True, help='The modularity of conservation')

    parser.add_argument('--genome', type=str, default='data/hg38.fasta', help='The Fasta file of genome')

    parser.add_argument('--gtf', type=str, default='data/Homo_sapiens.Ensembl.GRCh38.82.gtf',
                        help='The gtf annotation file. e.g., hg38.gtf')

    parser.add_argument('--bigwig', type=str, default='data/hg38.phyloP20way.bw',
                        help='conservation scores in bigWig file format')

    parser.add_argument('--positive_bed', type=str, default='data/circRNA_dataset.bed',
                        help='BED input file for circular RNAs for training, it should be like:chromosome    start    end    gene')

    parser.add_argument('--negative_bed', type=str, default='data/negative_dataset.bed',
                        help='BED input file for other long non coding RNAs for training, it should be like:chromosome    start    end    gene')
    parser.add_argument('--testing_bed', type=str, default='data/test.bed',
                        help='BED input file for testing data, it should be like:chromosome    start    end    gene')
CircDC.py  

3. Output explaining
The output file (in ".csv" format) can be found in results folder, which including sequence number, sequence_id, predicted probability and pedicted result.

