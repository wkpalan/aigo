from pylab import *

import sys, os

from AIGO import logger

from AIGO.ReferenceSet import RefSet
from AIGO.FunctionalAnnotation import FuncAnnot

from AIGO.Analyse import AnalyseFA
from AIGO.Plot    import PlotFA

from AIGO.go.OBO import readGOoboXML

from AIGO.utils.File import createDir
from AIGO.utils.Execute import batchExecute

from AIGO.Similarity import GOSet_Similarity, GOSet_PWSimilarity


def compareSimilarity(projectDir):
    """
    This function compare Similarity measures.
    """

    projectName="simPipeline"
    organism="bovine"

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

    allPipeName = ["AFFY", "B2G"]
    allFileType = allPipeName
    
    pipeline=dict()
    for pipeName, fileName, fileType in zip(allPipeName, allFileName, allFileType):
        FA=FuncAnnot(pipeName, refSet, G, organism=organism)
        FA.read(fileName, fileType=fileType)
        pipeline[pipeName]=FA

    #-----------------------------------------------
    
    #Analyse Functional annotations
    analyseFA = AnalyseFA()
    batchList=["removeUnconnected"]
    batchExecute(batchList, analyseFA, [pipeline[pipeName] for pipeName in allPipeName])


    #Compute information content
    logger.info("=================================================")
    logger.info("Computing Information Content")
    allIC=dict()
    for pipeName in allPipeName:
        FA=pipeline[pipeName]
        logger.info("\t%s" % FA.name)
        allIC[pipeName]=dict()
        for a in FA.G.aspect:
            allIC[pipeName][a]=dict()
            for go in FA.GOtoGP[a]:
                n=len(FA.GOtoGP[a][go])
                for ans in FA.G.ancestors(FA.G.get_intid(go)):
                    allIC[pipeName][a][ans]=allIC[pipeName][a].get(ans, 0) + n

        for a in FA.G.aspect:
            if len(allIC[pipeName][a].values())==0:
                continue
            m=max(allIC[pipeName][a].values())
            for go in allIC[pipeName][a]:
                allIC[pipeName][a][go]=-1. * log(1.*allIC[pipeName][a][go]/m)


    #Compare coherence of biological process annotation sets in AFFY given by three different similarity metrics
    logger.info("=================================================")
    logger.info("Computing functional coherence of biological process annotation sets in AFFY given by three different similarity metrics")
    aspect="biological_process"
    pipeName="B2G"
    FA=pipeline[pipeName]
    logger.info("\tGS2")
    allGS2=[mean(GOSet_Similarity(G, FA.GPtoGO[aspect][gp], metric="GS2")) for gp in FA.GPtoGO[aspect] if len(FA.GPtoGO[aspect][gp]) > 1]
    logger.info("\tCzekanowskiDice")
    allCD=[mean(GOSet_Similarity(G, FA.GPtoGO[aspect][gp], metric="CzekanowskiDice")) for gp in FA.GPtoGO[aspect] if len(FA.GPtoGO[aspect][gp]) > 1]
    logger.info("\tResnik")
    allResnik=[mean(GOSet_Similarity(G, FA.GPtoGO[aspect][gp], metric="Resnik", IC=allIC[FA.name])) for gp in FA.GPtoGO[aspect] if len(FA.GPtoGO[aspect][gp]) > 1]
    logger.info("\tCorrelation between GS2 and CzekanowskiDice : %.2f" % corrcoef(allGS2, allCD)[0][1])
    logger.info("\tCorrelation between CzekanowskiDice and Resnik: %.2f" % corrcoef(allCD, allResnik)[0][1])


    #Compare molecular function  annotation sets in AFFY and B2G using three similarity metrics
    logger.info("=================================================")
    logger.info("Comparing molecular function annotation sets in AFFY and B2G using three different similarity metrics")
    aspect="molecular_function"
    commonGene=set(pipeline["AFFY"].GPtoGO[aspect].keys()).intersection(pipeline["B2G"].GPtoGO[aspect].keys())
    logger.info("\tProcessing %d genes" % len(commonGene))
    allGS2, allCD, allResnik=list(), list(), list()
    for gp in commonGene:
        GO1=pipeline["AFFY"].GPtoGO[aspect][gp]
        GO2=pipeline["B2G"].GPtoGO[aspect][gp]
        
        allGS2.append(GOSet_PWSimilarity(G, GO1, GO2, metric="GS2")[0])        
        allCD.append(GOSet_PWSimilarity(G, GO1, GO2, metric="CzekanowskiDice")[0])
        allResnik.append(GOSet_PWSimilarity(G, GO1, GO2, metric="Resnik", IC=allIC[FA.name])[0])
    logger.info("\tCorrelation between GS2 and CzekanowskiDice : %.2f" % corrcoef(allGS2, allCD)[0][1])
    logger.info("\tCorrelation between CzekanowskiDice and Resnik: %.2f" % corrcoef(allCD, allResnik)[0][1])
    
if __name__ == '__main__':
    try:
        import psyco
        psyco.full()
    except ImportError:
        pass
    
    projectDir="../../example-data"

    compareSimilarity(projectDir)

