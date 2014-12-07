galpy-paper-figures
====================

:bold:Figures from the galpy paper (`Bovy 2015 <link>`__)

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
(figures up to and including 24, excluding code examples).

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

We can generate figures 15 and 16 using

.. code-block:: none

   python figure15+16.py figure15.png figure16.png

which prints the mean deviation in the radial and vertical
frequencies, and in the corresponding angles.

Similarly, figures 17 and 18 are produced by the following

.. code-block:: none

   python figure17+18.py figure17.png figure18.png

which prints the actions and the deviations in the actions,
frequencies, and angles. The calculation with
actionAngleIsochroneApprox takes a *very* long time (~10 hours). You
can create a version with coarser time sampling by editing the line
that says ``tts= ts[::1]`` to skip more values in ``ts`` (e.g., ``tts=
ts[::20]``, which only takes about half an hour).

Figure 19, which displays the focal length to use when using the
Staeckel approximation of the actions for MWPotential2014, can be
reproduced by

.. code-block:: none

   python figure19.py figure19.sav figure19.png

The savefile contains a pickle of the 2D array of focal lengths that
is displayed. See the ``figure19.py`` code for how to read this and
what the grid on which it is calculated is. The code prints the radius
of a circular orbit for each L grid point. This code also takes a long
time: about 2.5 hours.

The two panels of figure 20 can be obtained as

.. code-block:: none

   python figure20.py figure20a.ps figure20b.ps

Again, my version of matplotlib has some weird issues with plotting
the black points in the top panel for PNG output, which is why this
command is written to produce PS figures.

The two panels of figure 21 can be created using

.. code-block:: none

   python figure21.py figure21a.png figure21b.png

Figure 22 can in principle be produced by doing

.. code-block:: none

   python figure22.py figure22.png

but this will take a very long time, as all of the corrections
corresponding to different iterations have to be computed (it does not
take *forever*...).

Figure 24 is created by running

.. code-block:: none

   python figure24.py figure24.png

This can also take a very long time to run, especially if the
necessary DF corrections for the Dehnen and Shu DFs have not been
calculated before; progress in going through the various DFs is
printed. This code creates two files containing pickles of the
asymmetric drift for all models and the Oort constants. See the code
for more information.

The remaining figures require so much computation time to run that it
is not particularly interesting to exactly reproduce them. 