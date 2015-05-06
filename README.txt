In this folder the macros used for L1 upgrade tau trigger are saved and mantained.
First creation: 6 May 2015
Many programs previously realized will be progressively migrated here.

Structure:

** L1LegacyTrigger: programs to analyze Run I trigger (for comparison with Stage 2)
    - make rate evaluation
    - produce filtered taus (i.e. apply match hps-gen-L1)

** Stage2Trigger:  programs to analyze Stage 2 merge trees

** Common: programs to globally analyze data from both Stage2 and RunI triggers
    - Evaluate resolutions

** plotters: make plots by combining histograms from all other programs