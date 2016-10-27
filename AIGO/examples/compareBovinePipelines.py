# -*- coding: utf-8 -*-
from pylab import *

import sys, os

from AIGO import logger

from AIGO.ReferenceSet import RefSet
from AIGO.FunctionalAnnotation import FuncAnnot

from AIGO.Analyse import AnalyseFA
from AIGO.Plot    import PlotFA
from AIGO.Compare import CompareFA
from AIGO.Randomize import RandomizeFA
from AIGO.Report import ReportFA

from AIGO.go.OBO import readGOoboXML

from AIGO.utils.File import createDir
from AIGO.utils.Execute import batchExecute


def compareBovinePipelines(projectDir):
    """
    This function compare the properties of 3 functional annotations for a Bovine array.
    """
    
    projectName="bovinePipeline"
    organism="bovine"

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("This function compare the properties of 3 functional annotations for a Bovine array.")
    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")

    logger.info("name of the project : %s " % projectName)

    #Read bovine microarray probe set to define the set of gene products
    fileName= "%s/ReferenceSet/%s.fasta" %(projectDir, organism)
    refSet=RefSet(organism=organism, fileName=fileName, refType="Fasta")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)

    #Read Functional annotations
    allFileName=list()
    allFileName.append("%s/Annotation/Affy_%s.na31.annot.csv" % (projectDir, organism))
    allFileName.append("%s/Annotation/B2G_%s.annot" % (projectDir, organism))
    allFileName.append("%s/Annotation/AID_%s.txt" % (projectDir, organism))

    allPipeName = ["AFFY", "B2G", "AID"]
    allFileType = allPipeName
    
    pipeline=dict()
    for pipeName, fileName, fileType in zip(allPipeName, allFileName, allFileType):
        FA=FuncAnnot(pipeName, refSet, G, organism=organism)
        FA.read(fileName, fileType=fileType)
        pipeline[pipeName]=FA

    #-----------------------------------------------
    
    #Analyse Functional annotations
    analyseFA = AnalyseFA()
    #batchList=["obsolete", "unconnected", "removeUnconnected", "coverage",  "richness", "numberAnnot", "coherence", "redundancy", "compactness", "specificity", "informationContent"]
    batchList=["obsolete", "unconnected", "removeUnconnected", "coverage",  "richness", "numberAnnot", "redundancy", "specificity", "informationContent"]
    batchExecute(batchList, analyseFA, [pipeline[pipeName] for pipeName in allPipeName])
    
    #How big are the largest annotation sets ?
    analyseFA.largestSet([pipeline[pipeName] for pipeName in allPipeName])
    logger.info("The largest sets of annotations are :")
    for pipeName in allPipeName:
        FA=pipeline[pipeName]
        logger.info("\t%d for %s" % (FA['largestSet']['All_aspects_of_GO'], FA.name))        

    #Plot statistics of Functional annotations
    outDir="%s/Graph/%s" % (projectDir, organism)
    createDir(outDir)
    plotFA = PlotFA(xlabel="Annotation pipelines", outDir=outDir, name=projectName, organism=organism, ext="png")
    batchExecute(batchList, plotFA, [pipeline[pipeName] for pipeName in allPipeName], doGrid=True)

    #batchList=["coherenceHisto2D", "numberAnnotHisto2D"]
    batchList=["numberAnnotHisto2D"]
    batchExecute(batchList, plotFA, [pipeline[pipeName] for pipeName in allPipeName], doGrid=True, tit="")
    
    #-----------------------------------------------

    #Compare  Functional annotations
    compareFA = CompareFA()
    batchList=["venn", "funcSim"]
    batchExecute(batchList, compareFA, [pipeline[pipeName] for pipeName in allPipeName])

    #Plot statistics of the comparison between Functional annotations
    batchList=["venn", "funcSymSim"]
    batchExecute(batchList, plotFA, compareFA, [pipeline[pipeName] for pipeName in allPipeName], tit="")

    #-----------------------------------------------
    #Export statistics to Excel
    outDir="%s/Export/%s" % (projectDir, organism)
    createDir(outDir)

    #exportList=["unconnected", "coverage",  "richness", "numberAnnot",  "coherence",  "compactness", "specificity", "informationContent", "redundancy"]
    exportList=["unconnected", "coverage",  "richness", "numberAnnot",  "specificity", "informationContent", "redundancy"]
    reportFA = ReportFA(outDir=outDir, name=projectName, organism=organism)
    reportFA.printStatistics([pipeline[pipeName] for pipeName in allPipeName] ,exportList)
    reportFA.saveStatistics([pipeline[pipeName] for pipeName in allPipeName] ,exportList)

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("")


def compareBovinePipelinesAndMerge(projectDir):
    """
    This function compare the properties of 3 functional annotations for a Bovine array + the intersection of the three + the union of the three
    """

    projectName="compareBovinePipelineAndMerge"
    organism="bovine"

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("This function compare the properties of:")
    logger.info("\t - 3 functional annotations for a Bovine array")
    logger.info("\t - the intersection of the 3 functional annotations")
    logger.info("\t - the union of the 3 functional annotations")
    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")

    logger.info("name of the project : %s " % projectName)

    #Read bovine microarray probe set to define the set of gene products
    fileName= "%s/ReferenceSet/%s.fasta" %(projectDir, organism)
    refSet=RefSet(organism=organism, fileName=fileName, refType="Fasta")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)

    #Read Functional annotations
    allFileName=list()
    allFileName.append("%s/Annotation/Affy_%s.na31.annot.csv" % (projectDir, organism))
    allFileName.append("%s/Annotation/B2G_%s.annot" % (projectDir, organism))
    allFileName.append("%s/Annotation/AID_%s.txt" % (projectDir, organism))

    allPipeName = ["AFFY", "B2G", "AID"]
    allFileType = allPipeName
    
    pipeline=dict()
    for pipeName, fileName, fileType in zip(allPipeName, allFileName, allFileType):
        FA=FuncAnnot(pipeName, refSet, G, organism=organism)
        FA.read(fileName, fileType=fileType)
        pipeline[pipeName]=FA

    #Union of the three  FAs
    FA=FuncAnnot("Union", refSet, G, organism=organism)
    for pipeName in ["AFFY", "B2G","AID"]:
        FA.add(pipeline[pipeName])
    pipeline["Union"]=FA

    #Intersection of the three  FAs
    FA=FuncAnnot("Intersection", refSet, G, organism=organism)
    FA.add(pipeline["AFFY"])
    for pipeName in ["B2G","AID"]:
        FA.inter(pipeline[pipeName])
    pipeline["Intersection"]=FA

    allPipeName.extend(["Union", "Intersection"])
    #-----------------------------------------------
    
    #Analyse Functional annotations
    analyseFA = AnalyseFA()
    #batchList=["obsolete", "unconnected", "removeUnconnected", "coverage",  "richness", "numberAnnot", "coherence", "redundancy", "compactness", "specificity", "informationContent"]
    batchList=["obsolete", "unconnected", "removeUnconnected", "coverage",  "richness", "numberAnnot", "redundancy", "specificity", "informationContent"]
    batchExecute(batchList, analyseFA, [pipeline[pipeName] for pipeName in allPipeName])

    #How big are the largest annotation sets ?
    analyseFA.largestSet([pipeline[pipeName] for pipeName in allPipeName])
    logger.info("The largest sets of annotations are :")
    for pipeName in allPipeName:
        FA=pipeline[pipeName]
        logger.info("\t%d for %s" % (FA['largestSet']['All_aspects_of_GO'], FA.name))
        
    #Plot statistics of Functional annotations
    outDir="%s/Graph/%s" % (projectDir, organism)
    createDir(outDir)
    plotFA = PlotFA(xlabel="Annotation pipelines", outDir=outDir, name=projectName, organism=organism, ext="png")
    batchExecute(batchList, plotFA, [pipeline[pipeName] for pipeName in allPipeName], doGrid=True)

    #batchList=["coherenceHisto2D", "numberAnnotHisto2D"]
    batchList=["numberAnnotHisto2D"]
    batchExecute(batchList, plotFA, [pipeline[pipeName] for pipeName in allPipeName], doGrid=True, tit="")

    #-----------------------------------------------

    #Compare  Functional annotations
    compareFA = CompareFA()
    batchList=["venn", "funcSim"]
    batchExecute(batchList, compareFA, [pipeline[pipeName] for pipeName in allPipeName])

    #Plot statistics of the comparison between Functional annotations
    batchList=["venn", "funcSymSim"]
    batchExecute(batchList, plotFA, compareFA, [pipeline[pipeName] for pipeName in allPipeName], tit="")

    #-----------------------------------------------
    #Export statistics to Excel
    outDir="%s/Export/%s" % (projectDir, organism)
    createDir(outDir)

    #exportList=["unconnected", "coverage",  "richness", "numberAnnot",  "coherence",  "compactness", "specificity", "informationContent", "redundancy"]
    exportList=["unconnected", "coverage",  "richness", "numberAnnot",  "specificity", "informationContent", "redundancy"]
    reportFA = ReportFA(outDir=outDir, name=projectName, organism=organism)
    reportFA.printStatistics([pipeline[pipeName] for pipeName in allPipeName] ,exportList)
    reportFA.saveStatistics([pipeline[pipeName] for pipeName in allPipeName] ,exportList)

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("")

def compareRandomizePipelines(projectDir):
    """
    This function compare the properties of 3 randomized functional annotations for a Bovine array. 
    """

    projectName="randomizePipeline"
    organism="bovine"

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("This function compare the properties of 3 randomized functional annotations for a Bovine array.")
    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")

    logger.info("name of the project : %s " % projectName)

    #Read rice microarray target sequence to define the set of gene products
    fileName= "%s/ReferenceSet/%s.fasta" %(projectDir, organism)
    refSet=RefSet(organism=organism, fileName=fileName, refType="Fasta")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)

    #Read Functional annotations
    allFileName=list()
    allFileName.append("%s/Annotation/Affy_%s.na31.annot.csv" % (projectDir, organism))
    allFileName.append("%s/Annotation/B2G_%s.annot" % (projectDir, organism))
    allFileName.append("%s/Annotation/AID_%s.txt" % (projectDir, organism))


    allPipeName = ["AFFY", "B2G", "AID"]
    allFileType = allPipeName

    pipeline=dict()
    for pipeName, fileName, fileType in zip(allPipeName, allFileName, allFileType):
        FA=FuncAnnot(pipeName, refSet, G, organism=organism)
        FA.read(fileName, fileType=fileType)
        pipeline[pipeName]=FA
    

    # Randomize FA
    randomizeFA = RandomizeFA()

    #-----------------------------------------------
    # Shuffle functional annotation
    batchList=["shuffleAnnotation"]
    batchExecute(batchList, randomizeFA, [pipeline[pipeName] for pipeName in allPipeName])

    #Analyse Functional annotations
    analyseFA = AnalyseFA()
    batchList=["coherence", "redundancy"]
    batchExecute(batchList, analyseFA, [pipeline[pipeName] for pipeName in allPipeName])

    #Export statistics to Excel
    outDir="%s/Export/%s" % (projectDir, organism)
    createDir(outDir)
    exportList=["coherence", "redundancy"]
    report = ReportFA(name="Randomize shuffle", outDir=outDir, organism=organism)
    report.printStatistics([pipeline[pipeName] for pipeName in allPipeName] ,exportList)
    report.saveStatistics([pipeline[pipeName] for pipeName in allPipeName] ,exportList)

    #-----------------------------------------------
    # Resample functional annotation
    batchList=["sampleAnnotation"]
    batchExecute(batchList, randomizeFA, [pipeline[pipeName] for pipeName in allPipeName])

    #Analyse Functional annotations
    #batchList=["obsolete", "unconnected", "removeUnconnected", "coverage", "richness", "numberAnnot", "coherence", "redundancy", "compactness", "specificity", "informationContent"]
    batchList=["obsolete", "unconnected", "removeUnconnected", "coverage", "richness", "numberAnnot", "redundancy", "specificity", "informationContent"]
    batchExecute(batchList, analyseFA, [pipeline[pipeName] for pipeName in allPipeName])

    #Export statistics to Excel
    outDir="%s/Export/%s" % (projectDir, organism)
    createDir(outDir)
    #exportList=["coverage",  "numberAnnot",  "richness", "coherence",  "compactness", "specificity", "informationContent", "redundancy"]
    exportList=["coverage",  "numberAnnot",  "richness", "specificity", "informationContent", "redundancy"]
    report = ReportFA(name="Randomize sample", outDir=outDir, organism=organism)
    report.printStatistics([pipeline[pipeName] for pipeName in allPipeName] ,exportList)
    report.saveStatistics([pipeline[pipeName] for pipeName in allPipeName] ,exportList)

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("")


def compareBovineAndRandom(projectDir):
    """
    This function compare the properties of 3 functional annotations for a Bovine array + a randomize version of Affymetrix functional annotations
    """
    
    projectName="BovineAndRandom"
    organism="bovine"

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("This function compare the properties of 3 functional annotations for a Bovine array + a randomize version of Affymetrix functional annotations.")
    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")

    logger.info("name of the project : %s " % projectName)

    #Read bovine microarray probe set to define the set of gene products
    fileName= "%s/ReferenceSet/%s.fasta" %(projectDir, organism)
    refSet=RefSet(organism=organism, fileName=fileName, refType="Fasta")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)


    #Read Functional annotations
    allFileName=list()
    allFileName.append("%s/Annotation/Affy_%s.na31.annot.csv" % (projectDir, organism))
    allFileName.append("%s/Annotation/B2G_%s.annot" % (projectDir, organism))
    allFileName.append("%s//Annotation/AID_%s.txt" % (projectDir, organism))
    allFileName.append("%s/Annotation/Affy_%s.na31.annot.csv" % (projectDir, organism))


    allPipeName = ["AFFY", "B2G", "AID", "resample"]
    allFileType = ["AFFY", "B2G", "AID", "AFFY"]
    
    pipeline=dict()
    for pipeName, fileName, fileType in zip(allPipeName, allFileName, allFileType):
        FA=FuncAnnot(pipeName, refSet, G, organism=organism)
        FA.read(fileName, fileType=fileType)
        pipeline[pipeName]=FA


    # Randimize FA
    randomizeFA = RandomizeFA()
    analyseFA = AnalyseFA()

    #-----------------------------------------------
    # Shuffle functional annotation
    batchList=["sampleAnnotation"]
    batchExecute(batchList, randomizeFA, [pipeline[pipeName] for pipeName in [ "resample"] ])
    
    batchList=["coherence", "redundancy", "numberAnnot"]
    batchExecute(batchList, analyseFA, [pipeline[pipeName] for pipeName in allPipeName])

    #Plot statistics of Functional annotations
    outDir="%s/Graph/%s" % (projectDir, organism)
    createDir(outDir)
    plotFA = PlotFA(xlabel="Annotation pipelines", outDir=outDir, name="Resample", organism=organism, ext="png")
    batchExecute(batchList, plotFA, [pipeline[pipeName] for pipeName in allPipeName], doGrid=True)

    batchList=["coherenceHisto2D", "numberAnnotHisto2D"]
    batchExecute(batchList, plotFA, [pipeline[pipeName] for pipeName in allPipeName], doGrid=True, tit="")

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("")

def worseFunctionalSimilarity(projectDir):
    """
    This function identifies the ten most different annotation sets between Affymetrix and Blast2GO for a Bovine array
    """
    
    from AIGO.Similarity import GOSet_PWSimilarity
    from itertools import izip

    projectName="bovinePipeline"
    organism="bovine"

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("This function identifies the ten most different annotation sets between Affymetrix and Blast2GO for a Bovine array")
    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")

    logger.info("name of the project : %s " % projectName)

    #Read bovine microarray probe set to define the set of gene products
    fileName= "%s/ReferenceSet/%s.fasta" %(projectDir, organism)
    refSet=RefSet(organism=organism, fileName=fileName, refType="Fasta")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)

    fileName="%s/Annotation/Affy_%s.na31.annot.csv" % (projectDir, organism)
    FA1=FuncAnnot("AFFY", refSet, G, organism=organism)
    FA1.read(fileName, fileType="AFFY")

    fileName="%s/Annotation/B2G_%s.annot" % (projectDir, organism)
    FA2=FuncAnnot("B2G", refSet, G, organism=organism)
    FA2.read(fileName, fileType="B2G")

    #Analyse Functional annotations
    analyseFA = AnalyseFA()
    batchExecute(["removeUnconnected"], analyseFA, [FA1, FA2])

    outDir="%s/Graph/%s/WorseFuncSim" % (projectDir, organism)
    createDir(outDir)

    N=10
    logger.info("=================================================")
    logger.info("Plotting the %d most dissimilar annotation sets" % N)
    logger.info("directory : %s" % outDir)
    for aspect in G.aspect:

        commonGene=set(FA1.GPtoGO[aspect].keys()).intersection(FA2.GPtoGO[aspect].keys())

        logger.info("%s : processing %d annotation sets " % (aspect, len(commonGene)))
                
        allD1, allD2=list(), list()
        for i,g in enumerate(commonGene):
            sim,l=GOSet_PWSimilarity(G, FA1.GPtoGO[aspect][g], FA2.GPtoGO[aspect][g])
        
            allD1.append(l[0])
            allD2.append(l[1])

        allD=map(lambda D:((array(D[0])+array(D[1]))/2.), izip(allD1,allD2))

        idx=argsort(allD)

        for i in arange(0,N):
            gp=list(commonGene)[idx[i]]

            figName="%s/%s_annotation_%s_from_%s_%s.png" % (outDir, aspect, gp, FA1.name, FA2.name)

            ttl="%s annotations of %s from %s (green) and %s (red) : Functional similarity = %.2f" % (aspect.replace("_", " "), gp, FA1.name, FA2.name, allD[idx[i]])
            FA1.G.compare_InducedGraph(FA1.GPtoGO[aspect][gp],FA2.GPtoGO[aspect][gp], figName=figName, ttl=ttl)

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("")


def GOFrequencyBovinePipelines(projectDir):
    """
    This function plot the frequency of GO terms in three bovine functional annotation
    """

    projectName="bovinePipeline"
    organism="bovine"

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("This function plot the frequency of GO terms from 3 functional annotations  for a Bovine array")
    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")

    logger.info("name of the project : %s " % projectName)

    #Read bovine microarray probe set to define the set of gene products
    fileName= "%s/ReferenceSet/%s.fasta" %(projectDir, organism)
    refSet=RefSet(organism=organism, fileName=fileName, refType="Fasta")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)

    #Read Functional annotations
    allFileName=list()
    allFileName.append("%s/Annotation/Affy_%s.na31.annot.csv" % (projectDir, organism))
    allFileName.append("%s/Annotation/B2G_%s.annot" % (projectDir, organism))
    allFileName.append("%s/Annotation/AID_%s.txt" % (projectDir, organism))

    allPipeName = ["AFFY", "B2G", "AID"]
    allFileType = allPipeName
    
    pipeline=dict()
    for pipeName, fileName, fileType in zip(allPipeName, allFileName, allFileType):
        FA=FuncAnnot(pipeName, refSet, G, organism=organism)
        FA.read(fileName, fileType=fileType)
        pipeline[pipeName]=FA    

    #----------------------------------------------
    #Plot frequency of GO terms in a radial grah
    outDir="%s/Graph/%s" % (projectDir, organism)
    logger.info("=================================================")
    logger.info("Plotting frequency of GO terms")
    logger.info("directory : %s" % outDir)    
    for aspect in G.aspect:
        logger.info("%s : " % aspect)
        
        A=None
        for pipeName in allPipeName:
            l=array([log(1+len(pipeline[pipeName].GOtoGP[aspect].get(go, []))) for go in G.get_NodesfromAspect(aspect)])
            l=l/max(l)*256.
            l=[int(round(n)) for n in l]
            freq=dict([(n,c) for n,c in zip(G.get_NodesfromAspect(aspect), l)])

            figName="%s/Frequency_%s_%s.png" % (outDir, pipeline[pipeName].name, aspect)
            A=G.plot_FrequencyGraph(aspect, freq, figName=figName, ttl="", graphviz=A)

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("")

        
if __name__ == '__main__':
    try:
        import psyco
        psyco.full()
    except ImportError:
        pass
    
    projectDir="../../example-data"

    try:
        #compareBovinePipelines(projectDir)
        compareBovinePipelinesAndMerge(projectDir)
        #compareRandomizePipelines(projectDir)
        #compareBovineAndRandom(projectDir)
        worseFunctionalSimilarity(projectDir)
        #GOFrequencyBovinePipelines(projectDir)
        
    except Exception as e:
        print "AIGO is stopping because :"
        print "\t%s" % e.__doc__

    
