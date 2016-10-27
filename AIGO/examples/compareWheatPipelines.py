from pylab import *

import sys, os

from AIGO import allAspect

from AIGO.ReferenceSet import RefSet
from AIGO.FunctionalAnnotation import FuncAnnot

from AIGO.Analyse import AnalyseFA
from AIGO.Compare import CompareFA
from AIGO.Plot    import PlotFA
from AIGO.Randomize import RandomizeFA
from AIGO.Report import ReportFA

from AIGO.go.OBO import readGOoboXML

from AIGO.utils.File import createDir
from AIGO.utils.Execute import batchExecute

def AffyCCAnnotations(FA, organism):
    from AIGO.pyGS2 import plot_InducedGraph
        
    s=set()
    for gp in FA.GPtoGO['cellular_component']:
        s = s | FA.GPtoGO['cellular_component'][gp]
    
    outDir="%s/Graph/%s/CCAnnotations/" % (projectDir, organism)
    createDir(outDir)
        
    fileName="%s/annotation.png" % (outDir)
    plot_InducedGraph(list(s),G,fileName=fileName, ttl="Affymetrix wheat Cellular Component")

def createColorMap():
    import matplotlib
    
    cdict=    {'blue': ((0.0, 1.0, 1.0),
                        (0.000001, 1.0, 1.0),
                        (0.5, 0.0, 0.0),
                        (1.0, 0., 0.)),
               'green': ((0.0, 1.0, 1.0),
                         (0.000001, 0.0, 0.0),
                         (0.5, 1.0, 1.0),
                         (1.0, 0.0, 1.0)),
               'red': ((0.0, 1.0, 1.0),
                       (0.000001, 0.0, 0.0),
                       (0.5, 0.0, 0.0),
                       (1.0, 1.0, 1.0))}
    
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    #    pcolor(rand(10,10),cmap=my_cmap)
    #    colorbar()
    return my_cmap


def surface2D(X, Y, xname="", yname="", tit="", cblabel="", bins=20, interp=None, xlog=False, ylog=False):

    my_cmap=createColorMap()

    fig=figure()
    ax=fig.gca()
    if xlog:
        ax.set_yscale('log')
    if ylog:
        ax.set_xscale('log')
    
    #range=[ min(min(X),min(Y)), max(max(X),max(Y))]
    z,x,y=histogram2d(X, Y, bins=bins )
    xr=abs(x[1]-x[0])/2.
    yr=abs(y[1]-y[0])/2.
    #im=imshow(transpose(z), origin='lower', cmap=my_cmap,  extent=[x[0]-xr, x[-1]+xr, y[0]-yr, y[-1]+yr], aspect='auto', interpolation=interp)
    im=imshow(transpose(z), origin='lower', cmap=my_cmap,  extent=[x[0], x[-1], y[0], y[-1]], aspect='auto', interpolation=interp)

    P=[(a+b)/2 for a,b in zip(x[:-1], x[1:])]

    allMean, allY =list(), list()
    Y=array(Y)
    idx=digitize(X, x)-1
    for idbin in arange(len(x)-1):
        foo=Y[where(idx==idbin)]
        allY.append(foo)
        allMean.append(mean(foo))
        
    #ax.plot(P,allMean, c='red', lw=2)
    ax.plot((x[0:-1]+x[1:])/2., allMean, ls="None", color="red", marker=".", visible=False)
    draw()
    label=[l.get_position()[0] for l in ax.get_xticklabels()]

    xlim=ax.get_xlim()
    ylim=ax.get_ylim()
    
    
    d=ax.boxplot(allY, sym='', whis=0, positions=P, widths=abs(x[0]-x[1]))
    for l in d['medians']: l.set_linewidth(2)
    for l in d['boxes']: l.set_color('black')

    ax.xaxis.set_ticks(label)
    ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    #for l in ax.get_xticklabels(): l.set_rotation(90)    
    
    xlabel(xname)
    ylabel(yname)
    cb=colorbar()
    cb.set_label(cblabel)
    title(tit)
    
    return ax

def readCOEX_Filter(fileName):
    
    f=open(fileName, "rb")

    allProbe=f.next()[:-1].split(' ')
    MAT=dict()
    for row in f:
        row=row[:-1].split(' ')
        pA,pB=sort([allProbe[int(row[0])],allProbe[int(row[1])]])
        MAT.setdefault(pA, dict()).update({pB:float(row[2])});
    f.close()

    return MAT


def compareCoexpression(allFA):
    outDir="%s/Graph/%s/Coexpression" % (projectDir, organism)
    createDir(outDir)
    fileName="%s/correlations.tab_Filter_1000_0.001" % outDir
    MAT=readCOEX_Filter(fileName)

    for FA in allFA:
    
        allCOEX,allSim=dict(),dict()
        for aspect in allAspect:
            if aspect == "All_aspects_of_GO":
                continue

            allCOEX[aspect] = [MAT[pA][pB] for pA in MAT for pB in MAT[pA] if pA in FA.GPtoGO[aspect] and pB in FA.GPtoGO[aspect]]
            allSim[aspect] = [G.GS2([G.GOtoInt(FA.GPtoGO[aspect][pA]), G.GOtoInt(FA.GPtoGO[aspect][pB])])[0] for pA in MAT for pB in MAT[pA] if pA in FA.GPtoGO[aspect] and pB in FA.GPtoGO[aspect]]

            ax=surface2D(allCOEX[aspect],allSim[aspect], interp="bicubic", cblabel="Number of probeset", bins=50)
            #ax.scatter(X,Y, facecolor='none', edgecolor='black', alpha=0.1)
            ax.set_xlabel("Coexpression")
            ax.set_ylabel("Semantic Distance")

            ax.xaxis.grid()
            ax.yaxis.grid()

            title("GO %s in %s" % (aspect.replace("_", " "), FA.name))
            figName="%s/%s_COEXvsSIM_%s.png" % (outDir, FA.name, aspect)
            savefig(figName)

    from scipy import stats
    allKeys=[(pA,pB)for pA in MAT for pB in MAT[pA] if MAT[pA][pB] > 0.95]
    for FA in allFA:
        print "========================================================="
        print "FA: %s" % FA.name

        print "   GO Aspect \t Mean (Std) \t 50% [5%,95%] \t N"
        
        for aspect in allAspect:
            if aspect == "All_aspects_of_GO":
                continue
        
            Coex = [MAT[pA][pB] for pA,pB in allKeys if pA in FA.GPtoGO[aspect] and pB in FA.GPtoGO[aspect]]
            Sim = [G.GS2([G.GOtoInt(FA.GPtoGO[aspect][pA]), G.GOtoInt(FA.GPtoGO[aspect][pB])])[0] for pA,pB in allKeys if pA in FA.GPtoGO[aspect] and pB in FA.GPtoGO[aspect]]

            print "%s \t %.2f (%.2f) \t %.2f [%.2f,%.2f]\t %d" % (aspect, mean(Sim), std(Sim), stats.scoreatpercentile(Sim, 50), stats.scoreatpercentile(Sim, 5), stats.scoreatpercentile(Sim, 95), len(Sim))
                

def worseSim_COPSAandB2G(projectDir, organism, FA1, FA2):
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

            figName="%s/%s_annotation_%s_from_%s_%s.png" % (outDir, aspect, gp, FA1.name, FA2.name)
            #figName="%s/%s_annotation_%s_from_%s_%s.pdf" % (outDir, aspect, gp, FA1.name, FA2.name)

            ttl="%s annotations of %s from %s (green) and %s (red) : Functional similarity = %.2f" % (aspect.replace("_", " "), gp, FA1.name, FA2.name, allD[idx[i]])
            #ttl=""
            FA1.G.compare_InducedGraph(FA1.GPtoGO[aspect][gp],FA2.GPtoGO[aspect][gp], figName=figName, ttl=ttl)
    

#-------------------------------------------------------------
def compareWheatPipelines(projectDir):
    """
    Add some comments here
    """

    organism="wheat"

    #Read rice microarray target sequence to define the set of gene products
    fileName= "%s/ReferenceSet/%s.fasta" %(projectDir, organism)
    refSet=RefSet(organism=organism, fileName=fileName, refType="Fasta")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)

    projectName="MATT"

    #Read Functional annotations
    allFileName=list()
    allFileName.append("%s/Annotation/blast2goPaths_fin_aracyc_%s_unionBest.tab" % (projectDir, organism))
    allFileName.append("%s/Annotation/pfam2goPaths2_%s_unionBest.tab" % (projectDir, organism))
    allFileName.append("%s/Annotation/pfam2goPaths2_%s_unionBest___blast2goPaths_fin_aracyc_%s_unionBest_merged.tab" % (projectDir, organism, organism))
    allFileName.append("%s/Annotation/COPSA_%s.tab" % (projectDir, organism))
    allFileName.append("%s/Annotation/Affy_%s.annot.csv" % (projectDir, organism))
    allFileName.append("%s/Annotation/B2G_%s.annot" % (projectDir, organism))


    allPipeName = ["Blast", "Pfam", "Merge", "COPSA", "AFFY", "B2G"]
    allFileType = ["GP2GO", "GP2GO", "GP2GO", "GP2GO", "AFFY", "B2G"]

    pipeline=dict()
    for pipeName, fileName, fileType in zip(allPipeName, allFileName, allFileType):
        FA=FuncAnnot(pipeName, refSet, G, organism=organism)
        FA.read(fileName, fileType=fileType)
        pipeline[pipeName]=FA

    #Analyse Functional annotations
    analyseFA = AnalyseFA()
    batchList=["obsolete", "unconnected", "removeUnconnected", "coverage",  "richness", "numberAnnot", "coherence", "redundancy", "compactness", "specificity", "informationContent"]
    batchExecute(batchList, analyseFA, [pipeline[pipeName] for pipeName in allPipeName])

    #Plot statistics of Functional annotations
    outDir="%s/Graph/%s" % (projectDir, organism)
    createDir(outDir)
    plotFA = PlotFA(xlabel="Annotation pipelines", outDir=outDir, name=projectName, organism=organism)
    batchExecute(batchList, plotFA, [pipeline[pipeName] for pipeName in allPipeName], doGrid=True)

    batchList=["coherenceHisto2D", "numberAnnotHisto2D"]
    batchExecute(batchList, plotFA, [pipeline[pipeName] for pipeName in allPipeName], doGrid=True)

    #Compare  Functional annotations
    compareFA = CompareFA()
    batchList=["venn", "funcSim"]
    batchExecute(batchList, compareFA, [pipeline[pipeName] for pipeName in ["COPSA", "AFFY", "B2G"]])

    #Plot statistics of the comparison between Functional annotations
    batchList=["venn", "funcSymSim"]
    batchExecute(batchList, plotFA, compareFA, [pipeline[pipeName] for pipeName in ["COPSA", "AFFY", "B2G"]])

    #-----------------------------------------------
    #Export statistics to Excel
    outDir="%s/Export/%s" % (projectDir, organism)
    createDir(outDir)

    exportList=["unconnected", "coverage",  "numberAnnot",  "richness", "coherence",  "compactness", "specificity", "informationContent", "redundancy"]
    reportFA = ReportFA(outDir=outDir, name=projectName, organism=organism)
    reportFA.printStatistics([pipeline[pipeName] for pipeName in allPipeName] ,exportList)
    reportFA.saveStatistics([pipeline[pipeName] for pipeName in allPipeName] ,exportList)
#-------------------------------------------------------------


#-------------------------------------------------------------
def compare_COPSAandB2G(projectDir):
    """
    Add some comments here
    """

    organism="wheat"

    #Read rice microarray target sequence to define the set of gene products
    fileName= "%s/ReferenceSet/%s.fasta" %(projectDir, organism)
    refSet=RefSet(organism=organism, fileName=fileName, refType="Fasta")

    #Read GO ontoloy
    fileName= "%s/OBO/go_daily-termdb.obo-xml" %(projectDir)
    G= readGOoboXML(fileName, force=False)

    projectName="MATT"

    #Read Functional annotations
    allFileName=list()
    allFileName.append("%s/Annotation/COPSA_%s.tab" % (projectDir, organism))
    allFileName.append("%s/Annotation/B2G_%s.annot" % (projectDir, organism))


    allPipeName = [ "COPSA", "B2G"]
    allFileType = ["GP2GO", "B2G"]

    pipeline=dict()
    for pipeName, fileName, fileType in zip(allPipeName, allFileName, allFileType):
        FA=FuncAnnot(pipeName, refSet, G, organism=organism)
        FA.read(fileName, fileType=fileType)
        pipeline[pipeName]=FA


    #COPSA annotations but only for GPs that are also annotated by B2G
    FA=FuncAnnot("COPSAandB2G", refSet, G, organism=organism)
    FA.add(pipeline["COPSA"])
    #So remove the GP that are not in B2G
    for aspect in G.aspect:
        copsaOnly=set(FA.GPtoGO[aspect].keys()).difference(pipeline["B2G"].GPtoGO[aspect].keys())
        FA.removeGP(copsaOnly, myAspects=[aspect])
    pipeline[FA.name]=FA
    allPipeName.append(FA.name)

    #B2G annotations but only for GPs that are also annotated by COPSA
    FA=FuncAnnot("B2GandCOPSA", refSet, G, organism=organism)
    FA.add(pipeline["B2G"])
    #So remove the GP that are not in B2G
    for aspect in G.aspect:
        b2gOnly=set(FA.GPtoGO[aspect].keys()).difference(pipeline["COPSA"].GPtoGO[aspect].keys())
        FA.removeGP(b2gOnly, myAspects=[aspect])
    pipeline[FA.name]=FA
    allPipeName.append(FA.name)

    #COPSA annotations only
    FA=FuncAnnot("COPSAonly", refSet, G, organism=organism)
    FA.add(pipeline["COPSA"])
    #So remove the GP that are in B2G
    for aspect in G.aspect:
        b2g=pipeline["B2G"].GPtoGO[aspect].keys()
        FA.removeGP(b2g, myAspects=[aspect])
    pipeline[FA.name]=FA
    allPipeName.append(FA.name)

    #B2G annotations only
    FA=FuncAnnot("B2Gonly", refSet, G, organism=organism)
    FA.add(pipeline["B2G"])
    #So remove the GP that are in COPSA
    for aspect in G.aspect:
        copsa=pipeline["COPSA"].GPtoGO[aspect].keys()
        FA.removeGP(copsa, myAspects=[aspect])
    pipeline[FA.name]=FA
    allPipeName.append(FA.name)
    
    #Analyse Functional annotations
    analyseFA = AnalyseFA()
    batchList=["obsolete", "unconnected", "removeUnconnected", "coverage",  "richness", "numberAnnot", "coherence", "redundancy", "compactness", "specificity", "informationContent"]
    batchList=["removeUnconnected"]
    batchExecute(batchList, analyseFA, [pipeline[pipeName] for pipeName in allPipeName])


    #Plot statistics of Functional annotations
    outDir="%s/Graph/%s" % (projectDir, organism)
    createDir(outDir)
    plotFA = PlotFA(xlabel="Annotation pipelines", outDir=outDir, name=projectName, organism=organism)
    batchExecute(batchList, plotFA, [pipeline[pipeName] for pipeName in allPipeName], doGrid=True)


    compareCoexpression([pipeline[name] for name in ["COPSAandB2G", "B2GandCOPSA", "COPSAonly", "B2Gonly"]])

if __name__ == '__main__':

    try:
        import psyco
        psyco.full()
    except ImportError:
        pass
    
    projectDir="../example-data"
    #compareWheatPipelines(projectDir)




