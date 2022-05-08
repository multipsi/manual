# Installing the program
## Installing binaries

This is the simplest way to obtain MultiPsi.
It makes use of Conda, which is an open-source package and environment management system that runs on all major operating systems.

Binaries are currently available for
- MacOS
- Linux

In time, it will also be made available for Windows 

It is enough that you install the minimal installer for conda named miniconda that includes only conda, Python, the packages they depend on, and a small number of other useful packages, including pip, zlib and a few others.
Retrieve miniconda from the following website

> <https://docs.conda.io/en/latest/miniconda.html>

Install the version for 64 bit computers that comes with Python (>=3.6).

Start a conda terminal, or Anaconda Powershell as it is referred to on a Windows system. Conda supports multiple *environments*
and you start in the one named `base` as is typically indicated by the prompt.
To create a new environment named `mtpenv` and install MultiPsi, Jupyter notebook and k3d (and package dependencies such as VeloxChem, NumPy and SciPy) into it, you enter the following command line statement

```
$ conda create -n mtpenv multipsi k3d jupyterlab -c veloxchem -c conda-forge
```

You can list your conda environments

```
$ conda env list
```

The activated environment will be marked with an asterisk (the `base` environment to begin with) and you can activate your new environment with the command

```
$ conda activate mtpenv
```

as should be indicated by getting a modified prompt.

To use the orbital viewer, inside this newly created environment you need to run the commands:

```
$ conda activate mtpenv
```

Inside this newly created environment, you should now be ready to start a Jupyter notebook with the command

```
$ jupyter nbextension install --py --user k3d
$ jupyter nbextension enable --py --user k3d
```

You are now ready to open the notebook using the command:

```
$ jupyter-notebook
```

which should open in your default web browser. A notebook allows for interactive execution of Python code written into cells. You should now be able to import the MultiPsi and VeloxChem modules in a cell:

```
import veloxchem as vlx
import multipsi as mtp
```

and start calculations. See the [eChem](https://kthpanor.github.io/echem) book for a multitude of examples.

## Installing from source
