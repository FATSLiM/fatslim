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

Please do not use code from this branch as it may break at anytime and burn you bad! You have been warned.

------
Status
------

You can follow the progression `here <https://github.com/FATSLiM/fatslim/projects/1>`_


.. _FATSLiM: http://fatslim.github.io/
.. _license: https://github.com/FATSLiM/fatslim/blob/master/LICENSE
.. _`develop branch`: https://github.com/FATSLiM/fatslim/tree/develop
.. _`master branch`: https://github.com/FATSLiM/fatslim/tree/master
.. _`http://pythonhosted.org/fatslim`: http://pythonhosted.org/fatslim
.. _MDAnalysis: https://www.mdanalysis.org/
.. _`Jupyter`: https://jupyter.org/