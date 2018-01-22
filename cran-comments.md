# Version 0.1.2


#Resubmission
CRAN Comments
** running examples for arch 'i386' ... [17s] NOTE
Examples with CPU or elapsed time > 10s
        user system elapsed
RDperm 11.62      0   11.65
** running examples for arch 'x64' ... [23s] NOTE
Examples with CPU or elapsed time > 10s
        user system elapsed
RDperm 16.41   0.02   16.42

Fixed by adding simpler examples and "\dontrun" to the more complex ones

Retested it in x86_64-w64mingw32 (64-bit) with success,
uploaded also to \url{https://win-builder.r-project.org/upload.aspx}, checks for arch 'i386' and arch 'x64' are ok

#Resubmission2
Added references to the description file
Added url to the references

#Resubmission3
References were wrong they wanted the link next to the paper.

#First Submission
## Test environments
* local OS X install, R 3.4.1
* x86_64-redhat-linux-gnu (64-bit), R 3.4.0

## R CMD check results
There were no ERRORs or WARNINGs.