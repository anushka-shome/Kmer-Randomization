# Kmer-Randomization

The following libraries are imported in the program: sys, random, and collections. All of these packages should be included if you have Python installed.

To run the program, download the kmer.py file to your device. Open your terminal and navigate to the folder where the kmer.py file is stored. Run the program using the following command line format:
```
python kmer.py <file_name>.fasta <k> <s/e>
```
- file_name.fasta is the fasta file containing your input sequence
- k determines the length of the kmer
- s/e refers to the type of randomization. 's' will call the sampling function and 'e' will call the function that finds the euler path in a DeBruijn graph

If the parameters you enter are incorrect, the following will be printed out:
```
Usage: python kmer.py <file_name>.fasta <k> <s/e>
If using Euler: DNA sequences only work with k > 5 and protein sequences only work with k > 2
Otherwise, program will be stuck in an infinite loop
```

There are some limitations to my program. While the sampling randomization works for any sequence of any value of k > 0, the euler function does not. 
- If the sequence is a DNA sequence, only k > 5 works
- If the sequence is a protein sequence, only k > 3 works
My program does not handle the event that the user inputs invalid k values for the euler function. Instead, the program will be stuck in an infinite loop.

The program writes into a file called randomized_seq.fasta, which will be created in the same folder as kmer.py. If you run the program once and then run it again, the contents of randomized_seq.fasta will be overwritten.

If an invalid file is given in file_name.fasta, an error statement is printed out and execution stops. Other errors are not handled besides an invalid command line, as mentioned above.

The output of the program is a table containing the kmer count for both the original and the randomized sequence.

The Example Sequences folder contains nucleotide and protein sequences. The nucleotide sequence files have a name beginning with an 'N' and the protein sequence files have a name beginning with a 'P'.
