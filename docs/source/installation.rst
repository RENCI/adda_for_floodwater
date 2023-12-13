###########################
Installation / Virtual Env
###########################

* Clone this repo.  This locations is referred to as **PATHTO** below.

* Make a python virtual environment (venv) with the requirements.txt file, called **adda**.  The actual name of the venv does not matter, as long as the name is speficied correctly in the **data_assimilation.yaml** file (see below).

  * conda create --name adda --file requirements.txt

* Several required packages are not available through conda.  These need to be pip-installed:

  * conda activate adda

  * pip install -r pip.reqs.txt

* Add ADDA paths to the virtual environment via the **PYTHONPATH** variable.

::

    conda activate adda
    conda env config vars set PYTHONPATH=PATHTO/adda_for_floodwater:PATHTO/adda_for_floodwater/adda
