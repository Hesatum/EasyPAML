      seqfile = C:\Users\User\Desktop\EasyPAML\exemplos_teste\amostras\OG0076119.onlyexons.fas
     treefile = final-tree.txt
      outfile = OG0076119.onlyexons_M8_results.txt
   
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report
      seqtype = 1              * Data type
        ndata = 1              * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = 1              * Remove sites with ambiguity data?
		
        model = 0         * Models for ω varying across lineages
	  NSsites = 8          * Models for ω varying across sites
    CodonFreq = 2        * Codon frequencies
	  estFreq = 0              * Use observed freqs or estimate freqs by ML
        clock = 0              * Clock model
    fix_omega = 0         * Estimate or fix omega
        omega = 0.5        * Initial or fixed omega
