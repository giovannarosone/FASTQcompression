# FASTQcompression

## install

```sh
git clone --recursive https://github.com/giovannarosone/FASTQcompression.git 
cd FASTQcompression
make
```

## run

Given a string collection in reads.fastq:

```sh
python3 FASTQCompression.py dataset/reads.fastq -o results/output --all
```

```sh
Sending logging messages to file: results/output.log
=== gsufsort ===
external/gsufsort/gsufsort dataset/reads.fastq --bwt --qs -o results/output
Elapsed time: 0.1522
=== smooth-qs ===
src/fq_compression -e results/output.bwt -q results/output.bwt.qs -f dataset/reads.fastq -o results/output.fq
Elapsed time: 0.5366
=== header ===
sed -n 1~4p dataset/reads.fastq > results/output.h
Elapsed time: 0.0121
=== gsufsort ===
external/gsufsort/gsufsort results/output.fq --bwt --qs -o results/output.fq
Elapsed time: 0.1224
=== compression ===
7z a -mm=PPMd results/output.fq.7z results/output.fq
7z a -mm=PPMd results/output.h.7z results/output.h
7z a -mm=PPMd results/output.fq.bwt.7z results/output.fq.bwt
7z a -mm=PPMd results/output.fq.bwt.qs.7z results/output.fq.bwt.qs
Elapsed time: 0.3376
=== results ===
Original:	2.44 MB
Compressed:	0.91 MB
Ratio = 0.37
=== gzip ===
gzip -9 -k -f dataset/reads.fastq
Compressed:	0.59 MB
Ratio = 0.24
Elapsed time: 0.4503
=== 7z ===
7z a -mm=PPMd dataset/reads.fastq.7z dataset/reads.fastq
Compressed:	0.43 MB
Ratio = 0.18
Elapsed time: 0.1437
```


## FASTQcompression

Summing up what FASTQcompression does:

0. (Pre-processing) It takes in input a FASTQ file and build the ebwt string (with DNA bases) and the associated permuted qs string (with quality scores) by using the tool gsufsort.

1. It processes by means of fq_compression the ebwt string and the associated permuted qs string. Note that fq_compression takes in input the original FASTQ file only for recovering original headers.

2. fq_compression uses the code/library in https://github.com/nicolaprezza/bwt2lcp to set two binary vectors: 
LCP_minima and LCP_Threshold, where
    - LCP_Threshold[i]=1 iff LCP[i]>=k (where k is a threshold value, default 16),
    - LCP_minima=1 iff LCP[i] is a local minimum.

3. fq_compression detects and analyzes a cluster on the basis of LCP_minima and LCP_Threshold:
base symbols in any cluster may be changed according to the most frequent symbol in their own cluster, while quality score values in clusters are always smoothed according to four different strategies (described below).

4. fq_compression builds a new FASTQ file by inverting the original ebwt string through the LF mapping and by using the original FASTQ file to recover the headers.
The modified DNA symbols (stored in teh auxiliary vector BWT_MOD) replace their original symbol in the new FASTQ file, while the modifified QS replace the original QS. 

5. The new FASTQ file output of fq_compression is now compressed by using PPMd algorithm (7zip with option -mm=PPMd) according to both strategies below:
   1) the whole new FASTQ file is compressed.
   2) the new FASTQ file is first processed by gsufsort to build the new ebwt string and its associated permuted qs string. Then, the ebwt string and the permuted qs string are compressed by PPMd.
   
Note that the output file of fq_compression may have as headers just '@' (no titles). Indeed, if we ignore header information coming from the original FASTQ file, we can build a new FASTQ file just considering:
- the original ebwt string,
- the modified base symbols (stored in BWT_MOD),
- the modifified quality scores.

There are some parameter we can set to compile fq_compression in order to change the QS smoothing approach.

1) The parameter M is to choose one option among:
    - smoothing QS with MAX_QS (M=0)
    - smoothing QS with Avg_QS (M=1)
    - smoothing QS with default value (M=2)
    - smoothing QS with Mean_Err (M=3)
  
2) The parameter B is to use Illumina 8 level binning (B=1).

The default parameters are M=0 and B=0.

For example, to smooth QS by using Mean_Err+Illumina_8_level_binning, compile fq_compression by

```sh
cd src
make M=3 B=1
```
To run only fq_compression use the command:

```sh
cd ..
./src/fq_compression -e dataset/example.fq.ebwt -q dataset/example.fq.ebwt.qs -f dataset/example.fq -o result.fq
```
where
- example.fq.ebwt is the ebwt string,
- example.fq.ebwt.qs is the associated permuted qs string,
- example.fq is the original FASTQ file,
- result.fq is a new FASTQ file (output).




## References:

```sh
@thesis{thesisNiccoli,
  author       = {Enrico Niccoli}, 
  title        = {Compressione di sequenze di DNA in formato FASTQ},
  school       = {University of Pisa},
  note         = {Supervisors: Giovanna Rosone and Veronica Guerrini}
}
```

## Acknowledgement

Thanks to Nicola Prezza for providing part of the code and to Felipe Louza for helpful discussions.

---
<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html
