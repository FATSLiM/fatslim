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


Python Package Index (aka PyPI)
+++++++++++++++++++++++++++++++

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
Unfortunately, if you installed FATSLiM as a regular users, you may need to add it to your ``PATH``.

On Mac OS and GNU/Linux
~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    It is assumed that your shell is BASH_.

.. _BASH: https://en.wikipedia.org/wiki/Bash_(Unix_shell)

Simply add the following line to your `.bashrc`_ file:

.. code-block:: bash


.. _PATH: https://en.wikipedia.org/wiki/PATH_(variable)
.. _`.bashrc`: http://www.gnu.org/software/bash/manual/html_node/Bash-Startup-Files.html
