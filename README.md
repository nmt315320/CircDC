# CircDC
1. Description
 Circular RNAs (circRNAs) can regulate microRNA activity and are related to various diseases, such as cancer. Functional research on circRNAs is the focus of scientific research. Accurate identification of circRNAs is important for gaining insight into their functions.


We developed a novel framework, CircDC, for classifying circRNAs from other lncRNAs. CircDC uses four different feature encoding schemes and adopts a multilayer convolutional neural network and bidirectional long short-term memory network to learn high-order feature representation and make circRNA predictions. 

3. Availability

2.1 Datasets and source code are available at:
https://github.com/nmt315320/CircDC.git.

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
Changing working dir to DeepSoluE-master_source_code, and then running the following command:

python CircDC.py -i testing.fasta -o results.csv

-i: name of input_file in fasta format # folder “sequence” is the default file path of the input_file

-o name of output_file # folder “results” is the default file path for result save.

3. Output explaining
The output file (in ".csv" format) can be found in results folder, which including sequence number, sequence_id, predicted probability and pedicted result.

