from pylab import *

import sys, os

from AIGO import allAspect

from AIGO.ReferenceSet import RefSet
from AIGO.FunctionalAnnotation import FuncAnnot

from AIGO.Analyse import AnalyseFA
from AIGO.Compare import CompareFA
from AIGO.Plot import PlotFA
from AIGO.Randomize import RandomizeFA
from AIGO.Report import ReportFA

from AIGO.go.OBO import readGOoboXML

from AIGO.utils.File import createDir
from AIGO.utils.Execute import batchExecute

        #-------------------------------------------------------------
def compareRiceAffymetrixReleases(projectDir):
    """
    This function compare the properties of 10 release of Affymetrix annotations for a Rice array. 
    """

    projectName="Affymetrix"
    organism="rice"

    #Read rice microarray target sequence to define the set of gene products
    fileName= "%s/ReferenceSet/%s.fasta" %(projectDir, organism)
    refSet=RefSet(organism, fileName, refType="Fasta")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)

    #Read 11 release of Affymetrix Functional annotations
    release=arange(20,32)

    Affy=dict()
    for r in release:
        FA=FuncAnnot(str(r), refSet, G, organism=organism)
        fileName= "%s/Annotation/Affy_%s.na%d.annot.csv" % (projectDir, organism, r)
        FA.read(fileName, fileType="AFFY")
        Affy[r]=FA
    
    #Analyse Functional annotations
    analyseFA = AnalyseFA()
    batchList=["obsolete", "unconnected", "removeUnconnected", "coverage",  "richness"]
    batchExecute(batchList, analyseFA, [Affy[r] for r in release])

    #Plot statistics of Functional annotations
    outDir="%s/Graph/%s" % (projectDir, organism)
    createDir(outDir)
    plotFA = PlotFA(xlabel="Affymetrix Release number", outDir=outDir, name=projectName, organism=organism, ext="pdf")
    batchExecute(batchList, plotFA, [Affy[r] for r in release], doGrid=True, lloc="upper right")

    #Compare  release 20  and 31
    compareFA = CompareFA()
    batchList=["venn", "funcSim"]
    batchExecute(batchList, compareFA, [Affy[20], Affy[31]])

    #Plot statistics of the comparison 
    batchList=["venn", "funcSymSim"]
    batchExecute(batchList, plotFA, compareFA, [Affy[20], Affy[31]], doGrid=True, tit="")

    #Find the worse semantic similarity between 20 and 31
    worseFunctionalSimilarity(projectDir, organism, Affy[20], Affy[31])

    #Study the evolution of Glutamine Synthetase (GS)  annotations
    GSAnnotations(projectDir, organism, G, Affy)


def worseFunctionalSimilarity(projectDir, organism, FA1, FA2):
    from AIGO.Similarity import GOSet_PWSimilarity
    from itertools import izip

    for aspect in allAspect:

        if aspect=="All_aspects_of_GO":
            continue

        commonGene=set(FA1.GPtoGO[aspect].keys()).intersection(FA2.GPtoGO[aspect].keys())

        allD1, allD2=list(), list()
        for i,g in enumerate(commonGene):
            sim,l=GOSet_PWSimilarity(FA1.G, FA1.GPtoGO[aspect][g], FA2.GPtoGO[aspect][g])

            allD1.append(l[0])
            allD2.append(l[1])

        allD=map(lambda D:((array(D[0])+array(D[1]))/2.), izip(allD1,allD2))

        idx=argsort(allD)

        outDir="%s/Graph/%s/WorseFuncSim" % (projectDir, organism)
        createDir(outDir)

        for i in arange(0,10):
            gp=list(commonGene)[idx[i]]

            #figName="%s/%s_annotation_%s_from_%s_%s.png" % (outDir, aspect, gp, FA1.name, FA2.name)
            figName="%s/%s_annotation_%s_from_%s_%s.pdf" % (outDir, aspect, gp, FA1.name, FA2.name)

            #ttl="%s annotations of %s from %s (green) and %s (red) : Functional similarity = %.2f" % (aspect.replace("_", " "), gp, FA1.name, FA2.name, allD[idx[i]])
            ttl=""
            FA1.G.compare_InducedGraph(FA1.GPtoGO[aspect][gp],FA2.GPtoGO[aspect][gp], figName=figName, ttl=ttl)


def GSAnnotations(projectDir, organism, G, Affy):
    """
    This function study the evolution of the Affymetrix annotations for 4 GS probesets in Rice
    """

    GSProbes=['Os.7909.1.S1_at', 'Os.12728.1.S1_at', 'Os.7879.1.S1_at', 'Os.48875.1.S1_at']

    for gsp in GSProbes:
        print "============================================================================="
        print "                        Annotation for probe %s                              " % gsp
        print "============================================================================="
        

        d=dict()
        for r in Affy:
            print "------------Release %d---------------" % r
            for aspect in Affy[r].GPtoGO:
                if Affy[r].GPtoGO[aspect].has_key(gsp):
                    print "\t%s : " % aspect, Affy[r].GPtoGO[aspect][gsp]
                else:
                    print "\tNo annotations in %s" % aspect
                    
            #Get the list of all GO terms
            l=list()
            for aspect in Affy[r].GPtoGO:
                if Affy[r].GPtoGO[aspect].has_key(gsp):
                    l.extend(Affy[r].GPtoGO[aspect][gsp])

            if len(l) > 0:
                t=tuple(sort(l))
                if not d.has_key(t):
                    d[t]=r


        for t in d:
            outDir="%s/Graph/%s/GSAnnotations/%s" % (projectDir, organism, gsp)
            createDir(outDir)
        
            figName="%s/release%d.png" % (outDir, d[t])
            G.plot_InducedGraph(t,figName=figName, ttl="Affymetrix release %d" % d[t])



if __name__ == '__main__':

    try:
        import psyco
        psyco.full()
    except ImportError:
        pass
    
    projectDir="../../example-data/"
    compareRiceAffymetrixReleases(projectDir)
    
