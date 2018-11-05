# KLIC - Kernel Learning Integrative Clustering

## Installing 

In order to install this package, you will need to have _Rmosek_ installed. 

### Installing Rmosek [Linux]

* Download mosek from https://www.mosek.com/downloads/ to any folder (usually `/home/<username>/bin`). We will refer to it as `<my-mosek-directory>`.
* Unzip it:
```bash
cd <my-mosek-directory>
tar -xvf mosektoolslinux64x86.tar.bz2
```
* Add the following to your `~/.bashrc` file:

```bash
export PATH=$PATH:<my-mosek-directory>/mosek/8/tools/platform/linux64x86/bin
```

* If you are eligible apply for a free personal academic license at https://www.mosek.com/products/academic-licenses/ You will receive it (immediately) via email. 
* Put the license file in `<my-mosek-directory>/mosek`
* Proceed with the installation of _mosek_

```bash
export PKG_MOSEKHOME=<my-mosek-directory>/mosek/8/tools/platform/linux64x86
export PKG_MOSEKLIB=mosek64
```
* Now that you have _mosek_ installed, you can open R and install the _Rmosek_ package with 
```R
install.packages("Rmosek", type="source", INSTALL_opts="--no-multiarch", 
                  repos="http://download.mosek.com/R/8")
```

### Installing klic

In R, type
```R 
library(devtools)
install_github("acabassi/klic")
library(klic)
```
