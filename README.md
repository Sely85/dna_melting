dna_melting
===========

Welcome to dna_melting code!

dna_melting gives information about the sequence (length, GC content, molecular weigth) and compute melting temperature using:
 1. Wallance rule
 2. Salt adjusted method
 3. Khandelwal method
 4. Breslauer method (profile in bre_melting_curve.out)
 5. SantaLucia method (profile in san_melting_curve.out)
 6. Sugimoto method (profile in sug_melting_curve.out)
 7. Consensus method


BUILD (Linux)
-------------
g++ -lm dna_melting.cpp -o dna_melting


USAGE
-----
./dna_melting <inputfile> 
 
The inputfile should contain the following lines
- sequence (5'-->3') 
- salt concentration [M] (deal [Na+] = 0.05 M)
- total nucleotide strand concentration [M] (ideal concentration 5e-8M) 

For further information please check the manual ("dna_melting_manual.pdf").


EXAMPLE
-------
As an example, the melting temperature of a S1S2 sequence (GCGTCATACAGTGC), at [Na+]=0.05M with [DNA]=5e-8M, can be computed as follows:

      ./dna_melting S1S2.inp

The output provides information on the sequence (GC content, molecular weigth) and estimates of melting temperature, using different methods. 
The extimated curves of melting, computed using Breslauer, SantaLucia and Sugimoto methods, are also computed. The gnuplot file "plot_curve.gnu" can be used to generate a graph of such curves ("melting_curves.eps").  