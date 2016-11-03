# -*- coding: utf-8 -*-
from pylab import *

import sys, os

from AIGO import allAspect,logger

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


def compareEvidence(projectDir):
    """
    This function compare electronically infered and manually curated annotations to experimental annotations
    """

    projectName="EvidenceCode"
    organism="allSpecies"

    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")
    logger.info("This function compare electronically infered and manually curated annotations to experimental annotations.")
    logger.info("◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦◦")

    logger.info("name of the project : %s " % projectName)

    refSet=RefSet(organism=organism)

    allOrg=["Arabidopsis_thaliana","Drosophila_melanogaster","Mycobacterium_tuberculosis_ATCC_25618","Schizosaccharomyces_pombe","Bos_taurus","Escherichia_coli_ATCC_27325","Mycobacterium_tuberculosis_Oshkosh","Caenorhabditis_elegans","Escherichia_coli_MG1655","Oryza_sativa","Synechocystis_sp","Candida_albicans_SC5314","Gallus_gallus","Pseudomonas_fluorescens_Pf-5","Danio_rerio","Homo_sapiens","Rattus_norvegicus"]

    for refOrg in allOrg:
        #Define the set of gene products
        fileName= "%s/EvidenceCode/%s/two_experimental_evidence.goa" %(projectDir, refOrg)
        refSet.add(fileName=fileName, refType="GAF")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)

    #Read all annotations
    fileType="GAF"

    evidenceCodes=["EXP2", "IC", "TAS","ISS", "NAS","IEA"]

    allFA=dict()

    #-----------------------------------------------
    #Read Functional annotations obtained by experiments
    pipeName="EXP2"
    EXP2=FuncAnnot(pipeName, refSet, G, organism=organism)
    for refOrg in allOrg:
        fileName= "%s/EvidenceCode/%s/two_experimental_evidence.goa" %(projectDir, refOrg)
        FA=FuncAnnot(pipeName, refSet, G, organism=refOrg)
        FA.read(fileName, fileType=fileType)
        EXP2.add(FA)
    allFA[pipeName]=EXP2

    #-----------------------------------------------
    #Read Functional annotations obtained by human curation
    for pipeName in ["IC", "TAS","ISS", "NAS"]:
        EV=FuncAnnot(pipeName, refSet, G, organism=organism)
        for refOrg in allOrg:
            fileName= "%s/EvidenceCode/%s/%s.goa" %(projectDir, refOrg, pipeName)
            if not os.path.exists(fileName):
                continue
            FA=FuncAnnot(pipeName, refSet, G, organism=refOrg)
            FA.read(fileName, fileType=fileType)
            EV.add(FA)
        allFA[pipeName]=EV


    #Merge FAs Assigned by Human Curator
    FA=FuncAnnot("AHC", refSet, G, organism=organism)
    #for evidence in ["IC", "TAS", "ISS", "NAS"]:
    for evidence in ["IC", "ISS", "NAS"]:
        FA.add(allFA[evidence])
    allFA["AHC"]=FA

    #-----------------------------------------------
    #Read Functional annotations obtained without human curation
    for pipeName in ["IEA"]:
        EV=FuncAnnot(pipeName, refSet, G, organism=organism)
        for refOrg in allOrg:
            fileName= "%s/EvidenceCode/%s/%s.goa" %(projectDir, refOrg, pipeName)
            if not os.path.exists(fileName):
                continue
            FA=FuncAnnot(pipeName, refSet, G, organism=refOrg)
            FA.read(fileName, fileType=fileType)
            EV.add(FA)
        allFA[pipeName]=EV

    #-----------------------------------------------
    listFA=["EXP2", "AHC", "IEA"]

    #Analyse Functional annotations
    analyseFA = AnalyseFA()
    batchList=["obsolete", "unconnected", "removeUnconnected", "coverage",  "richness", "numberAnnot", "coherence", "redundancy", "removeRedundancy", "compactness", "specificity", "informationContent"]
    batchExecute(batchList, analyseFA, [allFA[evidence] for evidence in listFA])

    #Plot statistics of Functional annotations
    outDir="%s/Graph/%s" % (projectDir, organism)
    createDir(outDir)
    plotFA = PlotFA(xlabel="Evidence Codes", outDir=outDir, name=projectName, organism=organism)
    batchExecute(batchList, plotFA, [allFA[evidence] for evidence in listFA], doGrid=True)

    batchList=["coherenceHisto2D", "numberAnnotHisto2D"]
    batchExecute(batchList, plotFA, [allFA[evidence] for evidence in listFA], doGrid=True)

    #Compare  Functional annotations
    compareFA = CompareFA()
    batchList=["venn", "funcSim"]
    batchExecute(batchList, compareFA, [allFA[evidence] for evidence in listFA])
    batchList=["recall", "precision"]
    batchExecute(batchList, compareFA, [allFA[evidence] for evidence in listFA])

    #Plot statistics of the comparison between Functional annotations
    batchList=["venn", "funcSymSim"]
    batchExecute(batchList, plotFA, compareFA, [allFA[evidence] for evidence in listFA])
    batchList=["recall", "precision"]
    batchExecute(batchList, plotFA, compareFA, [allFA[evidence] for evidence in listFA])

    #-----------------------------------------------
    #Export statistics to Excel
    outDir="%s/Export/%s" % (projectDir, organism)
    createDir(outDir)

    exportList=["unconnected", "coverage",  "richness", "numberAnnot",  "coherence",  "compactness", "specificity", "informationContent", "redundancy"]
    reportFA = ReportFA(outDir=outDir, name=projectName, organism=organism)
    reportFA.printStatistics([allFA[evidence] for evidence in listFA] ,exportList)
    reportFA.saveStatistics([allFA[evidence] for evidence in listFA] ,exportList)

    #-----------------------------------------------
    # Invididual contributions of evidence codes
    logger.info("=================================================")
    logger.info("Invididual contributions of evidence codes")
    contribution=dict()
    for ec in ["IC", "TAS","ISS", "NAS", "IEA"]:
        contribution[ec]=set([(gp,go) for aspect in ["cellular_component", "molecular_function", "biological_process"] for gp in allFA[ec].GPtoGO[aspect] for go in allFA[ec].GPtoGO[aspect][gp] ])


    total_Annotation=sum([len(contribution[ec]) for ec in ["IC", "TAS","ISS", "NAS"] ])
    
    for ec in ["IC", "TAS","ISS", "NAS"]:
        logger.info("\t%.02f %% of the annotations are supported by %s" % (100. * len(contribution[ec]) / total_Annotation, ec))

    batchList=["recall", "precision"]
    batchExecute(batchList, compareFA, [allFA[evidence] for evidence in ["EXP2", "ISS", "TAS", "NAS", "IC", "AHC", "IEA"]])
    logger.info("Done")
    #-----------------------------------------------
    #Plotting precision and recall for each evidence code
    logger.info("=================================================")
    logger.info("Plotting precision and recall for each evidence code")
    reference="EXP2"
    #plotEvidence=["AHC", "IEA"]
    plotEvidence=["ISS", "TAS", "NAS", "IC", "AHC", "IEA"]
    
    evidenceMarker=dict(zip(plotEvidence, ['s', 'd', 'D', '*', 'p', 'h']))
    evidenceSize=dict(zip(plotEvidence, [8, 8, 8, 8, 15, 15]))
    aspectColor=dict(zip(allAspect, ["blue", "green", "red", "cyan"]))
    
    
    fig=figure(figsize=(8,8))
    for evidence in plotEvidence:

        for aspect in allAspect:
            if aspect == "All_aspects_of_GO":
                continue

            allX = compareFA['recall'][aspect][ (evidence, reference)].values()
            allY = compareFA['precision'][aspect][ (evidence, reference)].values()

            meanX = mean(allX)
            errX  = std(allX)/sqrt(len(allX))
            meanY = mean(allY)
            errY  = std(allY)/sqrt(len(allY))
            
            errorbar(meanX, meanY, xerr=errX, yerr=errY, alpha=0.9, hold=True, mfc=aspectColor[aspect], ecolor=aspectColor[aspect], marker=evidenceMarker[evidence], ms=evidenceSize[evidence])

    xlabel("Verspoor Hierarchical Recall")
    ylabel("Verspoor Hierarchical Precision")

    allMarker=['o', 'o', 'o', 's', 'd', 'D', '*', 'p', 'h']
    allColor=["green", "red", "cyan", "white", "white", "white", "white", "white", "white"]
    allLabel=[aspect.replace("_"," ") for aspect in allAspect if not aspect == "All_aspects_of_GO"]
    allLabel.extend(plotEvidence)

    
    foo=[ Line2D(arange(5), arange(5), ls='-', marker=m, color=c, label=l) for m,c,l in zip(allMarker, allColor, allLabel)]
    leg=legend(foo, allLabel, loc="upper left", numpoints=1)
    leg.legendPatch.set_alpha(0.5)

    grid()

    outDir="%s/Graph/%s" % (projectDir, organism)
    createDir(outDir)
    figName="%s/PrecisionVSRecall.png" % outDir
    savefig(figName)
    logger.info("Done")
    #-----------------------------------------------


        

if __name__ == '__main__':

    try:
        import psyco
        psyco.full()
    except ImportError:
        pass

    projectDir="."

    try:
        compareEvidence(projectDir)
        
    except Exception as e:
        print "AIGO is stopping because :"
        print "\t%s" % e
        print "\t%s" % e.__doc__
