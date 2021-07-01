net session>nul 2>&1

:main
set LSTC_LICENSE=local
set LSTC_FILE=C:\LSDYNA\program\lstc_file
"C:\ArupLocal\GitLab\Analytical_Models\SRA\Dyna_exe\ls-dyna_smp_d_R910_winx64_ifort131" i=Pushover_Bay15.key g=Pushover_Bay15.ptf f=Pushover_Bay15.thf o=Pushover_Bay15.otf ncpu=4
exit