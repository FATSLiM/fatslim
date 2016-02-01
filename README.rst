=======
FATSLiM
=======

`FATSLiM`_ means "**F**ast **A**nalysis **T**oolbox for **S**imulation of **Li**pid **M**embranes".
As you can imagine, the goal is to provide ~~the Ultimate Question of Life, the Universe, and Everything~~ a fast and efficient tool to extract membrane-related physical properties from Molecular Dynamics simulations of lipid bilayers.
To see what `FATSLiM`_ is able to do, please visit the homepage at http://fatslim.github.io/.

------------
Installation
------------

FATSLiM can be installed using pip via the following command:

::
    pip install fatslim


Alternatively, installation can be done directly from source code which is hosted in a git repository at https://github.com/FATSLiM/fatslim and is available under the GNU General Public License, version 3 (see `license`_).
You can then clone the repo:

::
    git clone https://github.com/FATSLiM/fatslim.git


And then compile and install using the usual ``setup.py``:

::
    cd fatslim/
    python setup.py install


Finally, it is always a good idea to run the ``self-test`` command to make sure that everything is OK:

::
    fatslim self-test


.. _FATSLiM: http://fatslim.github.io/
.. _license: https://github.com/FATSLiM/fatslim/blob/master/LICENSE
