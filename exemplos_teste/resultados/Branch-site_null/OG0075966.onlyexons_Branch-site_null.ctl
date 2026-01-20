      seqfile = C:\Users\User\Desktop\EasyPAML\exemplos_teste\amostras\OG0075966.onlyexons.fas
     treefile = labeled.nwk
      outfile = OG0075966.onlyexons_Branch-site_null_results.txt
   
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report
      seqtype = 1              * Data type
        ndata = 1              * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = 1              * Remove sites with ambiguity data?
		
        model = 2         * Models for ω varying across lineages
	  NSsites = 2          * Models for ω varying across sites
    CodonFreq = 7        * Codon frequencies
	  estFreq = 0              * Use observed freqs or estimate freqs by ML
        clock = 0              * Clock model
    fix_omega = 1         * Estimate or fix omega
        omega = 0.5        * Initial or fixed omega
