# FASTQcompression

## install

```sh
git clone --recursive https://github.com/giovannarosone/FASTQcompression.git 
cd FASTQcompression/external/gsufsort/
make DNA=1 TERMINATOR=0
cd ../../src
make
cd ../
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

Per compilarlo usare g++.
L'eseguibile si può lanciare usando il comando
```sh
./fastqcompression -e dataset/example.fq.ebwt -q dataset/example.fq.ebwt.qs -f dataset/example.fq -o result.fq
```

dove
example.fq.ebwt è la BWT
example.fq.ebwt.qs sono i QS permutati in accordo alla BWT
example.fq è il file originario
result.fq è il fastq di output.
