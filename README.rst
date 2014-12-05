galpy-paper-figures
====================

**Figures from the galpy paper (`Bovy 2015 <link>`__)**

AUTHOR
-------

Jo Bovy - bovy at ias dot edu

DEPENDENCIES
-------------

`galpy <https://github.com/jobovy/galpy>`__ and its dependencies.

PURPOSE
-------

This repository contains code to produce many of the figures in the
paper describing the `galpy <https://github.com/jobovy/galpy>`__ code
(figures up to and including 21, excluding code examples).

Figure 6 can be produced by calling (change the extension of the figure
filename to save as a different format)

.. code-block:: none
   python figure6.py figure6.png

Figures 7, 8, and 9 and the entries in Table 1 are produced by calling

.. code-block:: none
   python fitMWPotential2014.py fitResults.sav --rotcurve=figure8.png --kzcurve=figure7a.png --termcurve=figure7b.png --potname=figure9a.png --densname=figure9b.png --tablename=table1.txt

This routine fits a bulge+disk+halo potential to the data described in
section 3.5 in the galpy paper, producing MWPotential2014. Play around
with it and add your own data! The labels on some of the plots may be
somewhat cut off when using .png, but not when using .ps.

The two panels of figure 10 and figure 12 are produced by calling

.. code-block:: none
   python figure10+12.py figure10a.png figure10b.png figure12.ps

(my current version of matplotlib has a weird behavior in that it does not plot anything if I use figure12.png).

Figure 13 takes a while to compute (~20 min.) and is produced by
running

.. code-block:: none
   python figure13.py figure13.png

The code spits out the time the integration takes for each integrator
and the mean dE/E.

Figure 14 can be reproduced by

.. code-block:: none
   python figure14.py figure14.png

which again prints the time the integration takes for each integrator
and the mean of the Jacobian minus one.