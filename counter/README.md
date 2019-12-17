# CFIR Index

## Authors
Muhammad Haseeb and Fahad Saeed

# What do you need?
1. make
2. GCC v5.4.0 or later with C++11+ support
2. Windows: MinGW (x86_64-6.3.0-win32-seh-rt_v5-rev2)

# Pre CFIR-Index Steps
1. Digest the proteome database using Protein Digestion Simulator or OpenMS.
2. Remove the redundant peptide sequences in the digested database using DBToolkit.
3. Optional: The decoy database can be generated and appended to the target database using DBToolkit.
4. Convert the MS/MS data in (mzML/MS2) format using msconvert.exe. The MS/MS data are extracted and processed by MSToolkit.

# The Sample Driver Application
The sample application shows the software pipeline for successfully incorporating CFIR-Index for peptide search. The application firsts initializes Peptide and Ions Index respectively. Then the raw fragment-ion match score based query is conducted against CFIR-Index. The driver application creates an array called Matches (99999 hits x 1000 queries) which is filled with index number of candidate PSMs using CFIR-Querying algorithm. The candidate PSMs can be further processed or filtered as required by the application. However, the sample application just discards the results and records the execution time only. There is another version of sample application; that counts the PSMs discovered per query peptide. That application can be used to measure stats and can be requested from authors by email if required. 

# Configure CFIR-Index & Sample Application
1. Configure the CFIR-Index parameters in /slmtransform/include/config.h: 
*Make sure that #define/undef WINDOWS is correct according to your Host OS else there will be Seg faults.*
2. Add the mods information in the ./slm.cpp (Sample Application) using the following format:

## Format
    <some code> 
    
    /* Initialize the vModInfo */
    vModInfo.num_vars = #; // Types of modifications been specified (max 7)
    vModInfo.vmods_per_pep = #; // Number of modified AAs per peptide
    
    /* List of Mods Info */
    vModInfo.vmods[0].aa_per_peptide = #; // Number of modified AAs of this modification type
    vModInfo.vmods[0].modMass = # * SCALE; // Mass of the specified modification, for example: 79.97
    vModInfo.vmods[0].residues[0] = 'S'; // modified residues list (max 4 per mod type allowed), for example S
    vModInfo.vmods[0].residues[1] = 'T'; // T
    vModInfo.vmods[0].residues[2] = 'Y'; // Y
    
    vModInfo.vmods[0].aa_per_peptide = #; // Number of modified AAs of this modification type
    vModInfo.vmods[0].modMass = 15.997 * SCALE; // Mass of the specified modification, for example: 15.997
    vModInfo.vmods[0].residues[0] = 'M'; // modified residues list (max 4 per mod type allowed), for example M

# Building CFIR-Index Sample Application
1. Open the Git Bash shell (Windows) or normal Terminal (Ubuntu).
1. Navigate to CFIR-Index home directory /cfir/
2. Execute the following command: `make`

# Running CFIR-Index Sample Application
1. Navigate to /cfir
2. ./SLMTransform.exe

# Please Note:
1. We are integrating a full file based parameter parser for configuration settings to runtime.
2. Max digested peptide mass allowed: 5000Da.
3. The CFIR-Index only returns the shared-peak count scores not the XCorr scores, which can be used in any way e.g. post filtered based on precursor masses or sequence tags, formally scored for PSMs, FDRed.
4. Please contact us if you experience any trouble or bugs in the software. Thanks. :)

# Please cite our work
For queries or questions about CFIR-Index, please contact: {fsaeed, mhaseeb}@fiu.edu. Thank you.
