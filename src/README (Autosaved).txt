file.gpr => file-hits.txt
- read gpr file
** delete rows with flag = -100
** create a new index (id, name)
- calculate ratio
-- using bg
** using regression fit, bg ~ fg + sqrt(fg) + ln(fg)
** explore a little bit in R
- calculate mean and std (masked)
- calculate z-score
** using ratio
** using log(ratio)
- print gpr augmented with z-score up to threshold
- print top hits, one row per IOH

pool-hits.txt => clone-hits.txt
- read file-to-pool map
- read top hits
- find hits shared by row and column
- print as 12x12 table
