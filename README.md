#  Library Design script

## Files needed
- Library_Design.py  
- functions.py  
- input_parameters.json  
- List_RT.csv  
- Primer_univ.csv  


## Instructions


### Edit input_parameters.json

Change the settings for the design of your library, if necessary:
- `resolution`: integer = Size of the loci in nucleotides (ex: 20000)  
- `startLib`: integer =  Start genomic coordinate of the 1st locus (ex: 16000000)  
- `nbrProbeByLocus`: integer = Number of primary probes per locus (ex: 150)  
- `nbrLociTotal`: integer = Total number of loci (ex: 30)  
- `PrimerU`: string = Choice of the pair of universal primers 'primer1', 'primer2' until 'primer8' (ex: 'primer1')  
- `nbrBcd_RT_ByProbe`: integer = Number of the same barcode per primary probe (max=5)  
- `bcd_RT_File`: string = 'Barcodes.csv' or 'List_RT.csv'  
- `ChromosomeFolder`: string = Folder where the file containing the sequences homologous to the genomic DNA is located   
(ex: '/mnt/PALM_dataserv/DATA/Commun/genomes/dm6/OligoMiner/dm6_balanced')  
- `ChromosomeFile`: string = Name of the file containing the sequences homologous to the genomic DNA (ex: 'chr2L.bed')  

### Run script

For this, open a terminal at the script location. Run the following command to run the script:

```bash
$ python3 -m Library_Design.py 
```

## Documentation

The script performs the following steps to:

### Calculation of the coordinates for each locus  
The script will calculate the coordinates of each locus according to the starting coordinates of the 1st locus (`startLib`), the resolution (`resolution`) and the total number of barcodes chosen (`nbrLociTotal`).  

### Selection of primary probe sequences for each locus  
Then the script will search all the primary probe sequences according to each locus coordinate in the ChromosomeFile.  
If the number of primary probes found is greater than the number of desired primary probes (`nbrProbeByLocus`), the script will choose (randomly among all sequences) the number of desired primary probes. Conversely, if the number of primary probes found is smaller than the number of desired primary probes (`nbrProbeByLocus`), the script will leave all primary probes found.  


### Completion of primary sequences with barcodes/RT and universal primers  
To the sequence of the primary probes, the script will add the sequences of the barcodes/RT according to each locus, then the sequence of the universal forward and reverse primers.  

### Verification of the homogeneity of the size of the probes  
If the size of the primary probes varies by more than 10%, the script will perform a 3' completion with random nucleotides  

### Writing the different files in a result folder:  

- `1_Library_details.txt`: All primary probe sequences are grouped by barcode, with corresponding information (start-end-Bcr/RT). Each primary probe is subdivided into distinct parts in order to separate the different regions (primerU, Bcd/RT, sequence homologous to genomic DNA)  

- `2_Full_sequence_Only.txt`: File containing all the sequences of the raw primary probes (without any information). This is the file that is used for ordering the library.  

- `3_Library_Summary.csv`: Table summarizing all the information about each barcode (locus number, start, end, Bcd/RT, PU.fw, PU.rev, number of probes per locus...)   

- `4-OutputParameters.json`: JSON dictionary containing the parameters used to produce the library.  
