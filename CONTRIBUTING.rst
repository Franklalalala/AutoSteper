============
Contributing
============

Welcome to ``AutoSteper`` contributor's guide.

This document focuses on getting any potential contributor familiarized
with the development processes, but `other kinds of contributions`_ are also
appreciated.

If you are new to using git_ or have never collaborated in a project previously,
please have a look at `contribution-guide.org`_. Other resources are also
listed in the excellent `guide created by FreeCodeCamp`_ [#contrib1]_.

Please notice, all users and contributors are expected to be **open,
considerate, reasonable, and respectful**. When in doubt, `Python Software
Foundation's Code of Conduct`_ is a good reference in terms of behavior
guidelines.


Issue Reports
=============

If you experience bugs or general issues with ``AutoSteper``, please have a look
on the `issue tracker`_. If you don't see anything useful there, please feel
free to fire an issue report.

.. tip::
   Please don't forget to include the closed issues in your search.
   Sometimes a solution was already reported, and the problem is considered
   **solved**.

New issue reports should include information about your programming environment
(e.g., operating system, Python version) and steps to reproduce the problem.
Please try also to simplify the reproduction steps to a very minimal example
that still illustrates the problem you are facing. By removing other factors,
you help us to identify the root cause of the issue.


Documentation Improvements
==========================

You can help improve ``AutoSteper`` docs by making them more readable and coherent, or
by adding missing information and correcting mistakes.

``AutoSteper`` documentation uses Sphinx_ as its main documentation compiler.
This means that the docs are kept in the same repository as the project code, and
that any documentation update is done in the same way was a code contribution.


When working on documentation changes in your local machine, you can
compile them using |tox|_::

    tox -e docs

and use Python's built-in web server for a preview in your web browser
(``http://localhost:8000``)::

    python3 -m http.server --directory 'docs/_build/html'


Code Contributions
==================


Submit an issue
---------------

Before you work on any non-trivial code contribution it's best to first create
a report in the `issue tracker`_ to start a discussion on the subject.
This often provides additional considerations and avoids unnecessary work.

Setup environment
---------------------

Prepare a conda environment:

.. code:: 

   conda create -n AutoSteper python=3.8
   conda activate AutoSteper

Collect packages
--------------------

To install from the source code, the AutoSteper package:

.. code:: 

   git clone https://github.com/Franklalalala/AutoSteper
   cd AutoSteper
   pip install . -e

The FullereneDataParser package:

.. code:: 

   git clone https://github.com/XJTU-ICP/FullereneDataParser
   cd FullereneDataParser
   pip install . -e


.. warning::

   `FullereneDataParser <https://github.com/XJTU-ICP/FullereneDataParser>`__ contains part of C++ code, 
   to properly install, an advanced compiler version is required. Simply load the highest available version of 
   compiler will avoid most of the problems.


To compile the usenauty project, please follow instructions in
`usenauty <https://github.com/Franklalalala/usenauty>`__.


Submit your contribution
------------------------

AutoSteper implements `GitHub
Actions <https://github.com/features/actions>`__ to perform version
control. To trigger an effective response, commit messages needs to stay
in line with `Conventional Commit
messages <https://www.conventionalcommits.org/>`__.


Troubleshooting
---------------

Load the highest version of the aviable compiler in your enviroment will solve most of the problems. If not, `lower the version
requirement <https://github.com/Franklalalala/usenauty#note-for-compiler-version>`__
may help.



.. [#contrib1] Even though, these resources focus on open source projects and
   communities, the general ideas behind collaborating with other developers
   to collectively create software are general and can be applied to all sorts
   of environments, including private companies and proprietary code bases.


.. <-- strart -->
.. todo:: Please review and change the following definitions:

.. |the repository service| replace:: GitHub
.. |contribute button| replace:: "Create pull request"

.. _repository: https://github.com/Franklalalala/AutoSteper
.. _issue tracker: https://github.com/Franklalalala/AutoSteper/issues
.. <-- end -->


.. |virtualenv| replace:: ``virtualenv``
.. |pre-commit| replace:: ``pre-commit``
.. |tox| replace:: ``tox``


.. _black: https://pypi.org/project/black/
.. _CommonMark: https://commonmark.org/
.. _contribution-guide.org: https://www.contribution-guide.org/
.. _creating a PR: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request
.. _descriptive commit message: https://chris.beams.io/posts/git-commit
.. _docstrings: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
.. _first-contributions tutorial: https://github.com/firstcontributions/first-contributions
.. _flake8: https://flake8.pycqa.org/en/stable/
.. _git: https://git-scm.com
.. _GitHub's fork and pull request workflow: https://guides.github.com/activities/forking/
.. _guide created by FreeCodeCamp: https://github.com/FreeCodeCamp/how-to-contribute-to-open-source
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _MyST: https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html
.. _other kinds of contributions: https://opensource.guide/how-to-contribute
.. _pre-commit: https://pre-commit.com/
.. _PyPI: https://pypi.org/
.. _PyScaffold's contributor's guide: https://pyscaffold.org/en/stable/contributing.html
.. _Pytest can drop you: https://docs.pytest.org/en/stable/how-to/failures.html#using-python-library-pdb-with-pytest
.. _Python Software Foundation's Code of Conduct: https://www.python.org/psf/conduct/
.. _reStructuredText: https://www.sphinx-doc.org/en/master/usage/restructuredtext/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _tox: https://tox.wiki/en/stable/
.. _virtual environment: https://realpython.com/python-virtual-environments-a-primer/
.. _virtualenv: https://virtualenv.pypa.io/en/stable/

.. _GitHub web interface: https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files
.. _GitHub's code editor: https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files
