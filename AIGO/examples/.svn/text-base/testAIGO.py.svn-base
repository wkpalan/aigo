from pylab import *

from AIGO.FunctionalAnnotation import RefSet
from AIGO.FunctionalAnnotation import FuncAnnot

from AIGO.Analyse import AnalyseFA
from AIGO.PlotGraph import MultiPlotFA
from AIGO.Randomize import RandomizeFA
from AIGO.Report import ReportFA

from AIGO.IO import readGOoboXML

from AIGO.Utils import createDir
from AIGO.Utils import batchExecute

def testAIGO():
    from AIGO.examples.compareBovinePipelines import compareRiceAffymetrixReleases

def testAIGO():
    from compareRiceAffymetrixReleases import compareRiceAffymetrixReleases
    print ("Testing compareRiceAffymetrixReleases")
    compareRiceAffymetrixReleases("../../example-data");

    from compareBovinePipelines import compareBovinePipelines
    print ("Testing compareBovinePipelines")
    compareBovinePipelines("../../example-data")
    
    from compareBovinePipelines import compareRandomizePipelines
    print ("Testing compareRandomizePipelines")
    compareRandomizePipelines("../../example-data")

    from compareBovinePipelines import compareBovineAndRandom
    print ("Testing compareBovineAndRandom")
    compareBovineAndRandom("../../example-data")

if __name__ == '__main__':

    try:
        import psyco
        psyco.full()
    except ImportError:
        pass
    
    testAIGO()



    
