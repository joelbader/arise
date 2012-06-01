deconv
copyright (c) 2012 joel s. bader
joel.bader@jhu.edu

deconv process
- use argparse parse command-line arguments
-- put defaults into argparse help
- move to Dropbox/deconv
-- remove data files to make directory smaller
* make into a webserver trigger, check for new directories

control wells
- read a gpr file that has seth's flag = -100 indicating a control well
- generate a table with these ids
? we need separate tables for IgG vs. IgM controls
? not sure which is our reference

file.gpr => file-hits.txt
- read gpr file
-- delete rows with flag = -100
-- read controls from a file, based on id
-- delete rows with fg or bg == 0
-- create a new index (id, name)
- calculate ratio
-- using bg
** using regression fit, bg ~ fg + sqrt(fg) + ln(fg)
** explore a little bit in R
-- using norm based on second wavelength
- calculate mean and std (masked)
- calculate z-score
-- using ratio
-- using ratio(635 medians)/ratio(532 medians) to correct for GST
-- using log(ratio), with or without norm
- print gpr augmented with z-score up to threshold
-- create a results directory if it doesn't exist
- print top hits, one row per IOH

pool-hits.txt => clone-hits.txt
-- autogenerate a pool-to-file map from the filenames
-- probably good to include the directory in the filename
- read file-to-pool map
- read top hits
- find hits shared by row and column
- print as 12x12 table
