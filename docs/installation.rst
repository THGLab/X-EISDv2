============
Installation
============

X-EISDv2 uses only Python-based APIs for which we expect it to run
native on any system Python can run.

Please note that SPyCi-PDB is required to automate and perform back-calculations
if you do not provide back-calculated data to compare with experimental data for your
conformational ensembles. Details for SPyCi-PDB
can be found `here <https://github.com/julie-forman-kay-lab/SPyCi-PDB>`_.

Full installation instructures are highlighted below.

From Source
-----------

Clone from the official repository::

    git clone https://github.com/THGLab/X-EISDv2

Navigate to the new ``X-EISDv2`` folder::

    cd X-EISDv2

Run the following commands to install ``xeisd`` dependencies if
Anaconda is used as your Python package manager::

    conda env create -f requirements.yml
    conda activate xeisd

.. note::
    If you don't use Anaconda to manage your Python installations, you can use
    ``virtualenv`` and the ``requirements.txt`` file following the commands:

    | ``virtualenv spycienv --python=3.9``
    | ``source spycienv/bin/activate``
    | ``pip install -r requirements.txt``

    If you have difficulties installing ``spycipdb``, raise an issue in the
    main GitHub repository, and we will help you.

Install ``spycipdb`` in development mode in order for your installation to be
always up-to-date with the repository::

    python setup.py develop --no-deps

.. note::
    The above applies also if you used ``virtualenv`` instead of ``conda``.

**Remember** to active the ``xeisd`` environment every time you open a new
terminal window, from within the repository folder, choose yours::

    # Installation with Anaconda
    conda activate xeisd

    # Installation with virtualenv
    source xeisdenv/bin/activate

To update to the latest version, navigate to the repository folder, activate the
``xeisd`` python environment as described above, and run the commands::

    git pull

    # if you used anaconda to create the python environment, run:
    conda env update -f requirements.yml

    # if you used venv to create the python environment, run:
    pip install -r requirements.txt  --upgrade

    python setup.py develop --no-deps

Your installation will become up to date with the latest developments.


Installing SPyCi-PDB
--------------------

.. note::
    Please visit the SPyCi-PDB repository for more detailed installation instructions.
    Below highlights only the base installation of SPyCi-PDB.

    The ``xeisd`` conda environment should not be active. Deactivate using::
        
        conda deactivate
    
    if you're using ``virtualenv``, remain in the environment.

Clone from the official repository::

    git clone https://github.com/julie-forman-kay-lab/SPyCi-PDB

Navigate to the new ``SPyCi-PDB`` folder::

    cd SPyCi-PDB

Run the following commands to install ``spycipdb`` dependencies if
Anaconda is used as your Python package manager::

    conda env update --name xeisd --file requirements.yml --prune
    conda activate xeisd
    python setup.py develop --no-deps
    
Run the following commands to install ``spycipdb`` dependencies if
virtualenv was used to install SPyCi-PDB::

    pip install -r requirements.txt
    python setup.py develop --no-deps

Go back to the ``X-EISDv2`` directory and reinstall ``xeisd``::

    cd ..
    python setup.py develop --no-deps
