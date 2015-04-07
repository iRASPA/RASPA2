RASPA2
======

A general purpose classical simulation package that can be used for the
simulation of molecules in gases, fluids, zeolites, aluminosilicates,
metal-organic frameworks, carbon nanotubes and external fields.

Installation
============

Like any other software library, RASPA2 should be installed through a package
manager. This can be done through pip.

```
pip install RASPA2
```

If you want to stay on the cutting edge of RASPA2 development, you can also
install directly from this git repo with `pip install git+https://github.com/numat/RASPA2`.

If you are unfamiliar with package managers, read the
[For New Coders](https://github.com/numat/RASPA2/wiki/For-New-Coders)
section of the RASPA2 [wiki](https://github.com/numat/RASPA2/wiki).

Usage
=====

There are currently two ways to use RASPA: through configuration files in the
command line, or through Python functions.

###Command Line

This approach uses RASPA script files and a set of directories to load
structures, gases, and forcefield files. It outputs a set of folders and files
which contain data about the simulation run. This is the original designed use
of RASPA, and can be used in conjunction with shell scripts for small sets of
long-running simulations. To see the help, type:

```bash
simulate -h
```

To write and configure simulation input files, read
[Docs/raspa.pdf](https://github.com/numat/RASPA2/blob/master/Docs/raspa.pdf).

###Python

The previous approach is useful for long-running jobs, but can become
cumbersome when you want to manage a large number of simulations. In this use
case, the Python wrapper provides a streamlined workflow.

RASPA's full functionality can be accessed in a Python script, which enables
simulation runs to be "glued" together with auxiliary workflow steps. For
example, to run a set of simulations across a logarithmic range of pressures,
parse the outputs for uptake data, and plot the results, use:

```python
import RASPA2

# Set up
gas = "CO2"
pressures = [1e4 * 10**(0.1 * i) for i in range(21)]

# Run
results = [RASPA2.run(my_structure, gas, temperature=298, pressure=pressure)
           for pressure in pressures]

# Parse
uptakes = [r["Number of molecules"][gas]
            ["Average loading absolute [cm^3 (STP)/cm^3 framework]"][0]
           for r in results]

# Plot
plot(pressures, uptakes)
```

When used in conjunction with supercomputer interfaces such as the
[IPython Notebook](http://ipython.org/notebook.html) and chemical data handling
libraries such as [open babel](http://openbabel.org/wiki/Main_Page), it is
possible to completely automate job distribution, cif formatting, charge
estimation, and more in a single script.

For more examples, read the
[workflow wiki](https://github.com/numat/RASPA2/wiki/Workflow).

Accessing Configuration Files
=============================

In order to edit the RASPA forcefields, you need to access the associated .def
files. The current path of the .def files can be discovered by typing
`raspa-dir` into a terminal. In order to change your terminal directory to
this, you can chain commands:

```
cd `raspa-dir`
```

Edit the files here, and then re-run your simulations. If you make any
improvements to the existing force fields, molecules, or structures, you should
consider contributing them back to RASPA by issuing a
[pull request](https://help.github.com/articles/using-pull-requests/).
