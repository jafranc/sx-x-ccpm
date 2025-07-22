.. _command:

How it is launched
=======================

This tool is design to be used _via_ command line.

.. code:: bash
    :number-lines:

    /path/to/bin/ccpm --image /path/to/full.tiff -u 1,2,3 -o /path/to/output
    /path/to/bin/ccpm --image /path/to/output/isoVal_3_cc_1.tiff -n /path/to/output/isoVal_mapping.csv -x 5.1 -o /path/to/output

The first command will extract for each non-solid phase the :math:`n` images isolating :math:`n` connected components.
It will also output a csv table with the values and localization of phase transitions `isoVal_mapping.csv`.
For instance, points at the triple line will be flag :math:`1+2+3=6` there.

Eventually the `-x` or `--cutoff` is for input of the value from which to threshold for selecting the face elements. These elements are the
ones selected for integration of the Gauss-Bonnet :ref:`model` integration.
