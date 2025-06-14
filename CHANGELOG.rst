
Changelog
=========

v0.2.4 (2025-06-10)
------------------------------------------------------------

v0.2.3 (2025-05-13)
------------------------------------------------------------

* Minor bugfix for parsing of data files for ``optimize`` module

v0.2.2 (2025-02-25)
------------------------------------------------------------

* Add option to treat PRE back-calculated data as intensity ratios

v0.2.1 (2025-02-25)
------------------------------------------------------------

* Overhaul of normal-loglike function for better error checking
* Switch score and RMSE calculations to ignore NaN values

v0.2.0 (2025-01-15)
------------------------------------------------------------

* Update syntax for NOE and PRE experimental values from ``value`` to ``dist_value``

v0.1.5 (2024-11-19)
------------------------------------------------------------

* Update default chemical shift back-calculator errors to the ones referred to in 10.1039/C9SC06561J

v0.1.4 (2023-08-03)
------------------------------------------------------------

* Restore CLI capabilities and linting

v0.1.3 (2023-08-03)
------------------------------------------------------------

* Fix bug associated with optimize

v0.1.2 (2023-07-10)
------------------------------------------------------------

* Fix bug associated with old numpy version

v0.1.1 (2023-04-12)
------------------------------------------------------------

* Fix scoring and optimization in Rh module

v0.1.0 (2023-01-27)
------------------------------------------------------------

* Subsetting subclient to create optimized ensembles from indices results

v0.0.10 (2023-01-27)
------------------------------------------------------------

* Optimization subclient functionalities

v0.0.9 (2022-10-04)
------------------------------------------------------------

* CS scoring functions clean and test
* Integrate CS BC by SPyCi-PDB (UCBShift) if needed

v0.0.8 (2022-09-30)
------------------------------------------------------------

* SAXS scoring functions clean and test
* Integrate SAXS BC by SPyCi-PDB if needed

v0.0.7 (2022-09-28)
------------------------------------------------------------

* JC and smFRET scoring functions
* Integrate JC and smFRET BC by SPyCi-PDB if needed

v0.0.6 (2022-09-26)
------------------------------------------------------------

* Fix tests

v0.0.5 (2022-09-26)
------------------------------------------------------------

* Update requirements.txt for tests
* Lint and passing general tests

v0.0.4 (2022-09-23)
------------------------------------------------------------

* Score module (#4)
* Main logic for scoring functions in core/components directory
* Installation of IDPConformerGenerator on top of X-EISDv2
* Update README
* Installation of SPyCi-PDB on top of X-EISDv2
* ASCII art

v0.0.3 (2022-08-30)
------------------------------------------------------------

* Installation instructions (#3)

v0.0.2 (2022-08-30)
------------------------------------------------------------

* Core CLI and libraries (#2)
* Rename docs for README and setup.py

v0.0.1 (2022-08-29)
------------------------------------------------------------

* Housekeeping items (#1)
* Building based on python-project-skeleton
* Renaming and changing tests structures
