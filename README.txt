===========
   AIGO
===========

AIGO is a python library for the Analysis and 
the Inter-comparison of Gene Ontology functional annotations.
see (http://code.google.com/p/aigo).

Created by Michael Defoin-Platel on 21/02/2010.
Copyright (c) 2010. All rights reserved.

Typical usage could look like this::

    #!/usr/bin/env python

    from AIGO import logger

    from AIGO.ReferenceSet import RefSet
    from AIGO.FunctionalAnnotation import FuncAnnot
    from AIGO.go.OBO import readGOoboXML

    from AIGO.Analyse import AnalyseFA
    from AIGO.Report import ReportFA

    from AIGO.utils.Execute import batchExecute

    refSet = RefSet(organism="platypus", fileName="platypus.refSet", refType="Text")
    G = readGOoboXML("go_daily-termdb.obo-xml")
    FA = FuncAnnot("platypusProject", refSet, G, organism="platypus")
    FA.read("platypus.gaf", "GAF") 

    analyseFA = AnalyseFA()

    analyseFA.largestSet([FA])
    logger.info("Largest sets of annotations:")
    logger.info("\t%d for %s" % (FA['largestSet']['All_aspects_of_GO'], FA.name))

    batchList=["coverage",  "richness", "numberAnnot", "redundancy", "specificity", "informationContent"]
    batchExecute(batchList, analyseFA, [FA])

    reportFA = ReportFA(outDir=None, name="platypusProject", organism="platypus")
    reportFA.printStatistics([FA] ,batchList)

Tests
=====

Run testAIGO.py in the tests directory


Requirements
==============

Running AIGO on windows
-------------------------
* The 2.6.5 Python interpreter for Windows page http://www.python.org/download.
* GTK+ runtime (recommend bundle), PyGTK, PyCairo? and PyGObject http://www.pygtk.org/downloads.html
* BioPython? http://biopython.org/wiki/Download
* NumPy? http://sourceforge.net/projects/numpy/files/NumPy
* matplotlib http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.0
* xlwt http://pypi.python.org/pypi/xlwt 

Optional :
* wxPython http://www.wxpython.org/download.php#binaries
* psyco http://sourceforge.net/projects/psyco/files
* RPy http://sourceforge.net/projects/rpy/files


Contributors
============
* Michael Defoin-Platel
* Matthew Hindle
