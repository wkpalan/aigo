from pylab import *

import sys, os

import AIGO

from AIGO.FunctionalAnnotation import RefSet, FuncAnnot

from AIGO.Analyse import AnalyseFA
from AIGO.Plot import PlotFA
from AIGO.Randomize import RandomizeFA
from AIGO.Report import ReportFA

from AIGO.IO import readGOoboXML

from AIGO.utils.File import createDir
from AIGO.utils.Execute import batchExecute


#-------------------------------------------------------------
def testSCOP(projectDir):
    """
    """

    projectName="SCOP"
    organism="scop"

    #Read rice microarray target sequence to define the set of gene products
    fileName= "%s/ReferenceSet/%s.txt" %(projectDir, organism)
    refSet=RefSet(organism, fileName, refType="Text")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)

    fileName="%s/Annotation/%s.txt" % (projectDir, organism)
    
    pipeName=projectName
    fileType="SCOP"
    
    FA=FuncAnnot(pipeName, refSet, G, organism=organism)
    FA.read(fileName, fileType=fileType)
    
    #Analyse Functional annotations
    analyseFA = AnalyseFA()
    batchList=["unconnected", "removeUnconnected", "coverage",  "richness", "numberAnnot", "coherence", "redundancy", "compactness", "specificity", "informationContent"]
    batchExecute(batchList, analyseFA, [FA])

    #Plot statistics of Functional annotations
    outDir="%s/Graph/%s" % (projectDir, organism)
    createDir(outDir)
    plotFA = PlotFA(xlabel="", outDir=outDir, name=projectName, organism=organism)
    batchExecute(batchList, plotFA, [FA])

    #-----------------------------------------------
    #Export statistics to Excel
    outDir="%s/Export/%s" % (projectDir, organism)
    createDir(outDir)

    exportList=["unconnected", "coverage",  "numberAnnot",  "richness", "coherence",  "compactness", "specificity", "informationContent", "redundancy"]
    reportFA = ReportFA( outDir=outDir, name=projectName, organism=organism)
    reportFA.printStatistics([FA] ,exportList)
    reportFA.saveStatistics([FA] ,exportList)

#-------------------------------------------------------------
if __name__ == '__main__':
    import psyco
    psyco.full()
    projectDir="../../example-data"
    testSCOP(projectDir)




    
