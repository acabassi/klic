[![DOI](https://zenodo.org/badge/134583683.svg)](https://zenodo.org/badge/latestdoi/134583683)

# KLIC - Kernel Learning Integrative Clustering

## Installing 

In order to install this package, you will need to have _Rmosek_ installed. 

### Installing Rmosek [Unix]

* Download mosek from https://www.mosek.com/downloads/ to any folder (usually `/home/<username>/bin`). We will refer to it as `<my-mosek-directory>`.
* Unzip it:
```bash
cd <my-mosek-directory>
tar -xvf mosektoolslinux64x86.tar.bz2
```
* Add the following to your `~/.bashrc` file:

```bash
export PATH=$PATH:<my-mosek-directory>/mosek/<mosek-version>/tools/platform/<platform>/bin
```

`mosek-version` at the time of writing is `9.1`, whereas `<platform>` will be something like `linux64x86`, or `osx64x86`, or similar, depending on the machine.

* If you are eligible apply for a free personal academic license at https://www.mosek.com/products/academic-licenses/ You will receive it (immediately) via email. 
* Put the license file in `<my-mosek-directory>/mosek`
* Proceed with the installation of _mosek_

```bash
export PKG_MOSEKHOME=<my-mosek-directory>/mosek/<mosek-version>/tools/platform/<platform>
export PKG_MOSEKLIB=mosek64
```
* Now that you have _mosek_ installed, you can open R and install the _Rmosek_ package with 
```R
install.packages("Rmosek", type="source", INSTALL_opts="--no-multiarch", 
                  repos="http://download.mosek.com/R/<mosek-version>")
```

### Installing klic

In R, type
```R 
library(devtools)
install_github("acabassi/klic")
library(klic)
```

### Examples

- Clustering tumour samples: https://github.com/acabassi/klic-pancan12
- Clustering genes: https://github.com/acabassi/klic-transcriptional-modules
