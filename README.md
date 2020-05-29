# FASTQcompression

Per compilarlo usare g++.
L'eseguibile si può lanciare usando il comando
```sh
./fq_compression -e esempio.fq.ebwt -q esempio.fq.ebwt.qs -f esempio.fq -o risultato.fq
```

dove
esempio.fq.ebwt è la BWT
esempio.fq.ebwt.qs sono i QS permutati in accordo alla BWT
esempio.fq è il file originario
risultato.fq è il fastq di output.
