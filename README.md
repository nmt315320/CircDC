# CircDC
1. Description
CircRNA is the precondition for its industrial application and functional interpretation. However, the formation of inclusion bodies is still an inevitable roadblock in protein science and industry, where only nearly a quarter of proteins can be successfully expressed in soluble form.

CircDC is developed to predicts circRNA using a long-short-term memory (BLSTM) network with hybrid features composed of physicochemical patterns and distributed representation of amino acids. Comparison results showed that the proposed model achieved more accurate and balanced performance than existing tools.

2. Availability

2.1 Datasets and source code are available at:
 https://github.com/wangchao-malab/DeepSoluE/.

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

