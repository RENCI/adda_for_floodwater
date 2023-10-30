###########################
Installation / Virtual Env
###########################

* Clone this repo.  This locations is referred to as **PATHTO** below.
* Make a python virtual environment (venv) with the requirements.txt file, called **adda**.  The actual name of the venv does not matter, as long as the name is speficied correctly in the **data_assimilation.yaml** file (see below).
* conda create --name adda --file requirements.txt
* Several required packages are not available through conda.  These need to be pip-installed:
* conda activate adda
* pip install -r pip.reqs.txt
* Add ADDA paths to the virtual environment.  Either:
* Make a conda.pth in "envs/adda/lib/python3.8/site-packages/" that contains:

::

   **PATHTO**/adda_for_floodwater
   **PATHTO**/adda_for_floodwater/adda

* Note that the path to the venv's site-packages location may different than that above.  Modify as needed.
* or (probably better):
* Add the PYTHONPATH variable to the conda env. Activate the **adda** environment, add path:

::

    conda activate adda
    conda env config vars set PYTHONPATH=PATHTO/adda_for_floodwater:PATHTO/adda_for_floodwater/adda
