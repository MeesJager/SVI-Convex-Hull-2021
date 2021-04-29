
@echo off
:label1
set /p parallel="Do you want to compile the program to run in parallel? [Y/N]"
if /i "%parallel%" == "Y" (
    call gcc "main.c" -lm -fopenmp -D__USE_MINGW_ANSI_STDIO -o "HullAnalysis.exe"
) else (
    call gcc "main.c" -lm -lgomp -D__USE_MINGW_ANSI_STDIO -o "HullAnalysis.exe"
)
set /p continue="Do you want to end the compilation? [Y/N]"
if /i "%continue%" == "N" goto label1