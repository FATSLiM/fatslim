=======
FATSLiM
=======

`FATSLiM`_ means "**\ F**\ ast **\ A**\ nalysis **\ T**\ oolbox for **\ S**\ imulation of **\ Li**\ pid **\ M**\ embranes".
As you can imagine, the goal is to provide ~~the Ultimate Question of Life, the Universe, and Everything~~ a fast and efficient tool to extract membrane-related physical properties from Molecular Dynamics simulations of lipid bilayers.
To see what `FATSLiM`_ is able to do, please visit the homepage at http://fatslim.github.io/.

--------------
Important Note
--------------

`FATSLiM`_ is currently undergoing a major rewrite (read it is rewritten from scratch) in progress to overcome some intrinsic limitations.
The major goals of this rewrite are:

  - To use MDAnalysis`_ as a backend for trajectory and topology files I/O.
  - To provide an external API to use FATSLiM as a library so it can be use in e.g. `Jupyter`_ notebooks.
  - To ease code maintainability.

Please check out the `full-rewrite branch <https://github.com/FATSLiM/fatslim/tree/full-rewrite>`_ to follow the progression.


------
Status
------

`Master branch`_:
-----------------

.. image:: https://travis-ci.com/FATSLiM/fatslim.svg?branch=master
    :target: https://travis-ci.com/FATSLiM/fatslim

.. image:: https://coveralls.io/repos/github/FATSLiM/fatslim/badge.svg?branch=master
    :target: https://coveralls.io/github/FATSLiM/fatslim?branch=master

`Develop branch`_:
------------------

.. image:: https://travis-ci.com/FATSLiM/fatslim.svg?branch=develop
    :target: https://travis-ci.com/FATSLiM/fatslim

.. image:: https://coveralls.io/repos/github/FATSLiM/fatslim/badge.svg?branch=develop
    :target: https://coveralls.io/github/FATSLiM/fatslim?branch=develop


------------
Installation
------------

FATSLiM can be installed using pip via the following command:

.. code::

    pip install fatslim


Alternatively, installation can be done directly from source code which is hosted in a git repository at https://github.com/FATSLiM/fatslim and is available under the GNU General Public License, version 3 (see `license`_).
You can then clone the repo:

.. code::

    git clone https://github.com/FATSLiM/fatslim.git


And then compile and install using the usual ``setup.py``:

.. code::

    cd fatslim/
    python setup.py install


Finally, it is always a good idea to run the ``self-test`` command to make sure that everything is OK:

.. code::

    fatslim self-test


-------------
Documentation
-------------

A full documentation is available at: `http://pythonhosted.org/fatslim`_


.. _FATSLiM: http://fatslim.github.io/
.. _license: https://github.com/FATSLiM/fatslim/blob/master/LICENSE
.. _`develop branch`: https://github.com/FATSLiM/fatslim/tree/develop
.. _`master branch`: https://github.com/FATSLiM/fatslim/tree/master
.. _`http://pythonhosted.org/fatslim`: http://pythonhosted.org/fatslim
.. _MDAnalysis: https://www.mdanalysis.org/
.. _`Jupyter`: https://jupyter.org/