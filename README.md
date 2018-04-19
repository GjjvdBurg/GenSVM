GenSVM C Package
================

Paper: [GenSVM: A Generalized Multiclass Support Vector 
Machine](http://jmlr.org/papers/v17/14-526.html) by G.J.J. van den Burg and 
P.J.F. Groenen (*Journal of Machine Learning Research*, 2016).

GitHub: 
[https://github.com/GjjvdBurg/GenSVM](https://github.com/GjjvdBurg/GenSVM)

GenSVM is also available in these languages:

Language | URL
:-------:|:-------:
<a href="https://github.com/GjjvdBurg/PyGenSVM"><img src="https://www.python.org/static/community_logos/python-logo-master-v3-TM.png" height="75"/></a> | [https://github.com/GjjvdBurg/PyGenSVM](https://github.com/GjjvdBurg/PyGenSVM)
<a href="https://github.com/GjjvdBurg/RGenSVM"><img src="https://www.r-project.org/Rlogo.png" height="75"/></a> | [https://github.com/GjjvdBurg/RGenSVM](https://github.com/GjjvdBurg/RGenSVM)

Introduction
------------

GenSVM is a general multiclass support vector machine, which you can use for 
classification problems with multiple classes. Training GenSVM in 
cross-validation or grid search setups can be done efficiently due to the 
ability to use warm starts.  See the 
[paper](http://jmlr.org/papers/v17/14-526.html) for more information, and 
Usage below for how to use GenSVM.

The library has support for datasets in 
[MSVMpack](https://members.loria.fr/FLauer/files/MSVMpack/MSVMpack.html) and 
[LibSVM/SVMlight](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) format, and can 
take advantage of sparse datasets. There is also preliminary support for 
nonlinear GenSVM through kernels.

For documentation on how the library is implemented, see the [Doxygen 
documentation available here](https://gjjvdburg.github.io/GenSVM/). There are 
also many unit tests, which you can use to further understand how the library 
works. For the latest version of the library you can view the [test coverage 
report](https://gjjvdburg.github.io/GenSVM/cover) online.

This is the C library for GenSVM that contains two executables for using the 
method. A Python package for GenSVM is available 
[here](https://github.com/GjjvdBurg/PyGenSVM). An R package for GenSVM is 
planned.  If you are interested in this, please express your interest for the 
R package [here](https://github.com/GjjvdBurg/GenSVM/issues/2).

Usage
-----

First, download and compile the library. Minimal requirements for compilation 
are a working BLAS and LAPACK installation, which you can likely obtain from 
your package manager. It is however recommended to use ATLAS versions of these 
libraries, since this will give a significant increase in speed. If you choose 
not to use ATLAS, remove linking with ``-latlas`` in the ``LDFLAGS`` variable 
in the Makefile.

Then, compile the library with a simple:

    make

If you like to run the tests, use ``make test`` on the command line. 

After successful compilation, you will have two executables ``gensvm`` and 
``gensvm_grid``. Type:

    ./gensvm

To get an overview of the command line options to the executable (similar for 
``gensvm_grid``).

The ``gensvm`` executable can be used to train a GenSVM model on a dataset 
with a single hyperparameter configuration, whereas the ``gensvm_grid`` 
executable can be used to run a grid search on a dataset.

Here's an example of using the ``gensvm`` executable on a single dataset, with 
some custom parameters:

    ./gensvm -l 1e-5 -k 1.0 -p 1.5 data/iris.train

This fits the model with regularization parameter ``1e-5``, Huber hinge 
parameter ``1.0`` and lp norm parameter ``1.5``, and default settings 
otherwise. On my computer this yields a model with 18 support vectors in about 
0.1 seconds. The ``gensvm`` executable can also be used to get predictions for 
a test dataset, if it is supplied as final argument to the command. In this 
case, predictions will be printed to stdout, unless an output file is 
specified with the ``-o`` option.

The ``gensvm_grid`` executable can be used to run a grid search on a dataset.
The input to this executable is a file (called a grid file), which specifies 
the values of the parameters. See the ``training`` directory for examples and 
the documentation [here](https://gjjvdburg.github.io/GenSVM/) for more info on 
the file format.  One important thing to note is that when the ``repeats`` 
field has a positive value, a so-called "consistency check" will be performed 
after the grid search has finished. This is a robustness check on the best 
performing configurations, to find the best overall hyperparameter 
configuration with the best performance and smallest training time. In this 
robustness check warm-starts are not used, to ensure the observations are 
independent measurements of training time.

Here's an example of running ``gensvm_grid`` without repeats on the iris 
dataset:

    ./gensvm_grid training/iris_norepeats.training

On my computer this runs in about 8 seconds with 342 hyperparameter 
configurations. Alternatively, if consistency checks are desired we can run:

    ./gensvm_grid training/iris.training

which runs the same grid search but also does 5 consistency repeats for each 
of the configurations with the 5% best performance. Note that the performance 
is measured by cross-validated accuracy scores. This example runs in about 13 
seconds on my computer.

Reference
---------

If you use GenSVM in any of your projects, please cite the GenSVM paper 
available at 
[http://jmlr.org/papers/v17/14-526.html](http://jmlr.org/papers/v17/14-526.html). 
You can use the following BibTeX code:

    @article{JMLR:v17:14-526,
      author  = {Gerrit J.J. van den Burg and Patrick J.F. Groenen},
      title   = {{GenSVM}: A Generalized Multiclass Support Vector Machine},
      journal = {Journal of Machine Learning Research},
      year    = {2016},
      volume  = {17},
      number  = {225},
      pages   = {1-42},
      url     = {http://jmlr.org/papers/v17/14-526.html}
    }

License
-------

    Copyright 2016, G.J.J. van den Burg.

    GenSVM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GenSVM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GenSVM. If not, see <http://www.gnu.org/licenses/>.

    For more information please contact:

    G.J.J. van den Burg
    email: gertjanvandenburg@gmail.com
