# Software environment

The verification machine used R 4.5.1, Python 3.13.7, macOS on arm64, Pandoc 3.8, and pdfTeX from TeX Live 2024. `renv.lock` records the R package versions and their dependencies. Python code uses only the standard library.

The bootstrap worker uses R's fork-based parallel processing, which supports macOS and Linux. On Windows use one worker; multiple workers require a Unix-like environment. The full pipeline has not been verified on Windows.

To prepare a separate R library, first install `renv` (the lockfile was written with version 1.1.5). From the package directory run:

```sh
Rscript -e 'renv::restore(project=".", lockfile="renv.lock", library="/path/to/private-r-library", prompt=FALSE)'
export R_LIBS_USER=/path/to/private-r-library
```

Choose a library outside the release directory. Restoring packages may need internet access, a compiler, and system libraries for XML, curl/OpenSSL, fonts and text shaping, and linear algebra. R, Pandoc, and LaTeX are system dependencies and are not installed by `renv`.

Document rendering requires Pandoc and a LaTeX installation with the packages listed in the two manuscript wrappers. Use an existing TeX Live or TinyTeX installation. No proprietary software or remote computation service is required.

The numerical stage took approximately 32 minutes with eight CPU workers when reusing the manuscript's fitted models, excluding setup and rendering. Fewer workers preserve the same seed schedule but take longer. Full model-search runtime has not been measured for this package. Allow several GB of private disk space for copied inputs, models, figures, and saved analysis state.

The recorded R package versions were restored into a separate library using existing installed packages. A clean operating-system installation and other platforms have not been tested. See the main README for the scope of numerical verification.
