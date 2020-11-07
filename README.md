# FASTQcompression

## install

```sh
git clone --recursive https://github.com/giovannarosone/FASTQcompression.git 
make
```

## run

Given a string collection in reads.fastq:

```sh
python3 compress.py dataset/reads.fastq -o results/output --all
```

```sh
Sending logging messages to file: results/output.log
=== gsufsort ===
external/gsufsort/gsufsort dataset/reads.fastq --bwt --lcp --qs -o results/output
Elapsed time: 0.1522
=== smooth-qs ===
src/fastqcompression dataset/reads.fastq -e results/output.bwt -q results/output.bwt.qs -f dataset/reads.fastq -o results/output.fq
Elapsed time: 0.3771
=== header ===
sed -n 1~4p results/output.fq > results/output.h
Elapsed time: 0.0049
=== streams ===
sed -n 2~4p dataset/reads.fastq > results/output.dna
sed -n 4~4p dataset/reads.fastq > results/output.qs
Elapsed time: 0.0096
=== compression ===
7z a -mx9 -mmt12 results/output.h.7z results/output.h
7z a -mx9 -mmt12 results/output.dna.7z results/output.dna
7z a -mx9 -mmt12 results/output.qs.7z results/output.qs
Elapsed time: 0.6872
=== results ===
Original:	2.44 MB
Compressed:	0.43 MB
Ratio = 0.18
=== gzip ===
gzip -9 -k -f dataset/reads.fastq
Compressed:	0.59 MB
Ratio = 0.24
Elapsed time: 0.5498
=== 7z ===
7z a -mx9 -mmt12 dataset/reads.fastq.7z dataset/reads.fastq
Compressed:	0.49 MB
Ratio = 0.20
Elapsed time: 0.7864
```


## FASTQcompression
There are some parameter we can set to compile FASTQCompression in order to change the QS smoothing approach.

1) The parameter M is to choose one option among:
  - smoothing QS with MAX_QS (M=0)
  - smoothing QS with Avg_QS (M=1)
  - smoothing QS with default value (M=2)
  - smoothing QS with Mean_Err (M=3)
  
2) The parameter B is to use Illumina 8 level binning (B=1).

The default parameters are M=0 and B=0.

For example, to smooth QS by using Mean_Err+Illumina_8_level_binning, compile FASTQcompression by

```sh
cd src
make M=3 B=1
```
To run only FASTQcompression use the command:

```sh
cd ..
./src/fastqcompression -e dataset/example.fq.ebwt -q dataset/example.fq.ebwt.qs -f dataset/example.fq -o result.fq
```
where
example.fq.ebwt is the ebwt string,
example.fq.ebwt.qs is the associated permuted qs string,
example.fq is the original FASTQ file,
result.fq is a new FASTQ file (output).


Summing up what FASTQcompression does:

1. It takes in input the ebwt string and the associated permuted qs string, which have been obtained from a FASTQ file as pre-processing. 
Note that it takes in input the original FASTQ file only for recovering headers.

2. It uses the code/library in https://github.com/nicolaprezza/bwt2lcp to set two binary vectors: 
LCP_minima and LCP_Threshold, where
  - LCP_Threshold[i]=1 iff LCP[i]>=k (where k is a threshold value, default 16),
  - LCP_minima=1 iff LCP[i] is a local minimum.

3. It detects and analyzes a cluster on the basis of LCP_minima and LCP_Threshold:
eBWT symbols in any cluster may be changed according to the most frequent symbol in their own cluster, while QS values in clusters are always smoothed according to the choosen strategy.

4. It builds a new FASTQ file by inverting the original ebwt string through the LF mapping and by using the original FASTQ file to recover the headers.
The modified DNA symbols (stored in BWT_MOD) replace their original symbol in the new FASTQ file, while the modifified QS replace the original QS.

5. It outputs the new FASTQ file.

The output file may have as headers just '@' (no titles). Indeed, if we ignore header information coming from the original FASTQ file, we can build a FASTQ file starting from
- the original ebwt string,
- the modified DNA symbols (stored in BWT_MOD),
- the modifified QS.
