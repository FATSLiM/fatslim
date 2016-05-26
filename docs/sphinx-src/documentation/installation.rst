Installing FATSLiM
==================

Requirements
------------

FATSLiM is written in Python and, consequently, a **Python environment with header files** is required. FATSLiM is
developed and daily tested with Python **2.7 to 3.5** on GNU/Linux, but it also works on Windows and Mac OS
(far less tested on both, though).

FATSLiM only depends on `NumPy`_, but, if you want to run the tests (recommended), `pytest`_ is also needed.

Computationally intensive tasks are written in `Cython`_. Cython is not a requirement *per se*
as C-translated code is provided but the latest version of Cython is recommended.

Obviously, a **C compiler** is required to compile C extensions. Additionally, openMP_ is strongly
recommended as it is needed to enable parallelism.

.. _NumPy: http://www.numpy.org/
.. _pytest: http://pytest.org/
.. _Cython: http://cython.org/
.. _openMP: http://openmp.org/

Pre-installation
----------------

This section describes briefly how to install all the requirements needed by FATSLiM.
One caveat though: the goal here is not to provide a complete "How to install Python and stuff on your OS" guide (which would be out of the scope of this documentation)
but rather to give rough guidelines.

GNU/Linux
+++++++++

Python is most probably already installed on your system (at least for the vast majority of `Linux distributions`_).
To install the rest of the requirements you should use the `package manager`_ shipped with your system.
Here will be described only a few examples that may need to be adapted for your specific distro.
In such case, please refer to the distribution's documentation and/or package list.

.. _Linux distributions: https://en.wikipedia.org/wiki/Linux_distribution
.. _package manager: https://en.wikipedia.org/wiki/Package_manager

.. note::

    In the following, the package manager commands are given as if they are executed as root_.
    Depending on your distro, you may also use sudo_ to run the commands (very probably available).
    Please refer to your distribution manual.

.. _root: https://en.wikipedia.org/wiki/Superuser
.. _sudo: https://en.wikipedia.org/wiki/Sudo

.. note::

    Please also note that the provided command install the default packages. Depending on your distro,
    the packages may correspond to Python version 2 or 3.

APT-based distributions (Debian, Ubuntu, Linux Mint)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You should use apt_ to install the requirements:

.. _apt: https://en.wikipedia.org/wiki/Advanced_Packaging_Tool

.. code-block:: bash

    apt-get install gcc python-dev python-numpy cython python-pytest

YUM-based distributions (Fedora, RHEL, CentOS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You should use yum_ to install the requirements:

.. _yum: https://en.wikipedia.org/wiki/Yellowdog_Updater,_Modified

.. code-block:: bash

    yum install gcc python-devel numpy Cython pytest

ZYpp-based distributions (OpenSUSE, SUSE)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You should use zypper_ to install the requirements:

.. _zypper: https://en.wikipedia.org/wiki/ZYpp

.. code-block:: bash

    zypper install gcc python-devel python-numpy python-Cython python-pytest

Pacman-based distributions (Arch Linux, Manjaro)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You should use pacman_ to install the requirements:

.. _pacman: https://en.wikipedia.org/wiki/Arch_Linux#Pacman

.. code-block:: bash

    pacman -S gcc python-numpy cython python-pytest

Mac OS X
++++++++

The latest version of Mac OS X (10.11 - El Capitan) comes with Python 2.7 pre-installed.
Yet, the shipped version may be outdated and it is generally recommended by "Mac OS pythonistas" to use
homebrew_ to install a fully-functional and up-to-date Python version.
A complete guide to do so is available here:

`http://docs.python-guide.org/en/latest/starting/install/osx/`_.

.. _homebrew: http://brew.sh/
.. _`http://docs.python-guide.org/en/latest/starting/install/osx/`: http://docs.python-guide.org/en/latest/starting/install/osx/

Windows
+++++++

Installing Python and python packages on Windows can be cumbersome. Thankfully, several `Python distributions`_ exist and contain Python plus widely-used packages such as Numpy_ or Matplotlib_.
Anaconda_ is such Python distribution which has the advantage to be oriented toward science and is then well-suited to run FATSLiM.

This distribution can be downloaded from `https://www.continuum.io/downloads`_.

.. _Python distributions: https://wiki.python.org/moin/PythonDistributions
.. _Anaconda: https://www.continuum.io/why-anaconda
.. _Matplotlib: http://matplotlib.org/
.. _`https://www.continuum.io/downloads`: https://www.continuum.io/downloads

Once Anaconda is installed, you can access your Python environment via a dedicated command prompt
called ``Anaconda prompt`` accessible from the Windows Start menu.

Installation
------------

The easiest, and then recommended, way to install FATSLiM is probably to use pip_ and the `Python Package Index`_ (aka PyPI).
Alternately, you can also install FATSLiM from GitHub_.

.. _pip: http://www.pip-installer.org/en/latest/index.html
.. _Python Package Index: https://pypi.python.org/pypi
.. _GitHub: https://github.com

.. note::

    Installation of FATSLiM inside a `Python virtual environment`_ is of course possible but it will not be described here as setting up such environment is out of this documentation's scope.
    If you want to take this path, just follow the system-wide installation instructions **once** your virtual environment is activated.

.. _`Python virtual environment`: http://docs.python-guide.org/en/latest/dev/virtualenvs/


Python Package Index (recommended)
++++++++++++++++++++++++++++++++++

Installation-ready packages for FATSLiM are hosted on PyPI_. You can then use pip_ to install FATSLiM, independently of your OS.

.. _PyPI: https://pypi.python.org/pypi/fatslim

User-specific installation
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you do not have a `administrator privileges`_ on your machine (or simply do not want to mess with your system),
you can install (or upgrade) FATSLiM in your home directory with the following command:

.. code-block:: bash

    pip install --user --upgrade fatslim

.. _administrator privileges: https://en.wikipedia.org/wiki/Superuser

System-wide installation
~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    This also works for an installation inside a virtual environment (without the need for admin rights)

If you want to install/upgrade FATSLiM for all users, run the following command as `root/superuser/administrator`_:

.. code-block:: bash

    pip install --upgrade fatslim

.. _root/superuser/administrator: https://en.wikipedia.org/wiki/Superuser


GitHub
++++++

The main source code repository for FATSLiM is hosted on GitHub_ at the following address:
`https://github.com/FATSLiM/fatslim`_.
This repository can be used as an alternative to PyPI.

.. _`https://github.com/FATSLiM/fatslim`: https://github.com/FATSLiM/fatslim

End user installation
~~~~~~~~~~~~~~~~~~~~~

For each FATSLiM release, the corresponding package is accessible via the `"releases"`_ page where you can download the source code as a compressed archive (.zip or .tar.gz).
Once you have downloaded and extracted the archive, go to the source folder (``fatslim-x.y.z``).

You may then install FATSLiM in user site-package (recommended):

.. code-block:: bash

    python setup.py install --user

**Or**\ , you can install FATSLiM for all users (this requires `administrator privileges`_):

.. code-block:: bash

    python setup.py install

.. _"releases": https://github.com/FATSLiM/fatslim/releases


Coder/advanced user installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

    This section is targeted to people with a bit of practice with git_.

If you are a code/git_ guru, you can clone the repository to update FATSLiM continously instead of waiting for new releases (But you probably know that!)
For that matter, you need to clone the `main repository`_ using, for instance, the following:

.. _git: https://git-scm.com/
.. _main repository: https://github.com/FATSLiM/fatslim

.. code-block:: bash

    git clone https://github.com/FATSLiM/fatslim.git


The repository contains two branches that you may use:

* ``master``: Same code as the last public release
* ``develop``: Code not officially released but validated by `unit testing`_

.. _unit testing: https://en.wikipedia.org/wiki/Unit_testing

You can use this local repository to run FATSLiM but you will need to run ``python setup.py build_ext -i`` every time you switch between branches or pull the code.
You may then run FATSLiM by running ``fatslim`` (just besides ``setup.py``).


Post-installation
-----------------

Adding FATSLiM to your PATH
+++++++++++++++++++++++++++

If you performed a system-wide installation, FATSLiM executable is most likely already in your PATH_ and the ``fatslim`` command should be available.
Unfortunately, if you installed FATSLiM as a regular user and the ``fatslim`` command is not available, you need to add it to your ``PATH``.

On Mac OS and GNU/Linux
~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    It is assumed that your shell is BASH_.

.. _BASH: https://en.wikipedia.org/wiki/Bash_(Unix_shell)

Simply add the following lines to your `.bashrc`_ file:

.. code-block:: bash

    # Add Python user directory to PATH
    PYTHON_USER_PATH=`python -c "import site; print site.getuserbase()"`/bin
    if [ -d "$PYTHON_USER_PATH" ]
    then
      export PATH=$PATH:$PYTHON_USER_PATH
    fi

On Windows
~~~~~~~~~~

If you installed FATSLiM using Anaconda_ and the ``pip install --upgrade``, the ``fatslim`` script
is located inside the folder ``%CONDA_DEFAULT_ENV%\Scripts`` which is already in your ``%PATH%``.
Unfortunately, Windows does not recognize shebang_ lines and then the script is not executable.
The trick is to create a `Batch file`_ that will actually launch FATSLiM. To do so, follow these steps:

.. _shebang: https://en.wikipedia.org/wiki/Shebang_(Unix)
.. _Batch file: https://en.wikipedia.org/wiki/Batch_file

1. Open a Anaconda prompt (from Start menu) and run to get the default folder for Anaconda:

.. code-block:: bat

    echo %CONDA_DEFAULT_ENV%

2. With the file explorer, go to the ``Scripts`` folder inside the Anaconda folder. You should see the ``fatslim`` script there.
   Create a file named ``fatslim.cmd`` and paste the following code inside this new file:

.. code-block:: bat

    @echo off
    echo About to run FATSLiM with the following arguments: %*
    pause
    python %CONDA_DEFAULT_ENV%\Scripts\fatslim %*

3. Verify that can run FATSLiM from the Anaconda prompt by running ``fatslim --version``. The output should be similar to the following:

.. code-block:: bat

    [Anaconda2] C:\Users\IEUser>fatslim --version
    About to run FATSLiM with arguments: --version
    Press any key to continue . . .
    FATSLiM - Fast Analysis Toolbox for Simulations of Lipid Membranes
    version 0.1.2
    Copyright (c) 2013-2016 Sebastien Buchoux <sebastien.buchoux@gmail.com>

    FATSLiM - Fast Analysis Toolbox for Simulations of Lipid Membranes

    FATSLiM version: 0.1.2
    Python version: 2.7.11 (C:\Users\IEUser\Anaconda2\python.exe)
    Cython version (file generation): 0.23.4
    Python compiler: MSC v.1500 64 bit (AMD64)
    CPU architecture: 64bit
    OpenMP: 8 CPUs (default number of threads: 8 - max: 8)
    NumPy version: 1.11.0
    'version' command executed in 31.000 ms (CPU)
    Goodbye!

    [Anaconda2] C:\Users\IEUser>


Enable autocompletion
+++++++++++++++++++++

On Mac OS and GNU/Linux
~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    Currently, only Bash shell is supported for autocompletion.
    (This is the default shell for Mac OS and the vast majority of Linux distros)

Autocompletion for FATSLiM is enabled by the ``fatslim-autocompletion.bash`` script which is shipped with FATSLiM
First thing is to locate the ``fatslim-autocompletion.bash`` which is coherent with the FATSLiM version. To do so, run:

.. code-block:: bash

    which fatslim

This will output the actual path to the ``fatslim`` executable. This path is formatted like ``[ROOT]/bin/fatslim`` where ``[ROOT]`` is the root directory where FATSLiM-related files are installed.
For example, if ``which fatslim`` returns ``/usr/local/bin/python``, it means that the root directory is ``/usr/local/``.
Once you have, you know the root directory you can enable autocompletion by running the command (of course you need to replace ``[ROOT]`` by the actual value):

.. code-block:: bash

    source [ROOT]/share/fatslim/fatslim-autocompletion.bash

This will enable autocompletion only for the current session: if you close the terminal, autocompletion will be disabled. To enable it by default, simply add the whole ``source`` command to your `.bashrc`_.

On Windows
~~~~~~~~~~

Unfortunately for Windows users, autocompletion with a command prompt is rather limited and à la bash fully-fledged autocompletion is not feasable on Windows (at least easily).




.. _PATH: https://en.wikipedia.org/wiki/PATH_(variable)
.. _`.bashrc`: http://www.gnu.org/software/bash/manual/html_node/Bash-Startup-Files.html
