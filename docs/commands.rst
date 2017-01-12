Commands
========

The Makefile contains the central entry points for common tasks related to this project.


 * `make help` shows available commands.

Environment
^^^^^^^^^^^
Before running the project, make sure all requirements are installed in a
virtual enviroment.


 * `make requirements` will install Python dependencies of the project.
 * `make update-requirements` will update all python dependencies to the
   latest compatible version.

Running
^^^^^^^

 * `make data` will download and construct all dataset.


Tools
^^^^^
* `make lint` will use flake8 to check code.
* `make clean` will remove all files that can be build.
