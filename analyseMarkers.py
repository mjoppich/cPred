from collections import defaultdict, Counter
import math
import sys, os
import argparse
import gzip
from natsort import natsorted
import urllib.request
import lxml.html

if __name__ == '__main__':

    requiredColumns = ["gene", "clusterID", "avg_logFC", "p_val_adj", "mean", "num", "anum"]

    parser = argparse.ArgumentParser(description='Name Seurat clusters')
    parser.add_argument('-i', '--markers', type=argparse.FileType('r'), required=True, help='output from Seurat FindMarkers or FindAllMarkers')
    parser.add_argument('-up', '--update-panglao', action='store_true', default=False, help='update panglao db file')
    parser.add_argument('-uc', '--update-cellmarkerdb', action='store_true', default=False, help='update panglao db file')

    parser.add_argument('-g', '--gene', default="gene", type=str, help="column containing cluster id value")
    parser.add_argument('-c', '--cluster', default="clusterID", type=str, help="column containing cluster id value")

    parser.add_argument('-l', '--logfc', default="avg_logFC", type=str, help="column containing cluster id value")
    parser.add_argument('-p', '--pvaladj', default="p_val_adj", type=str, help="column containing cluster id value")
    parser.add_argument('-e', '--expr-mean', default="mean", type=str, help="column containing cluster id value")

    parser.add_argument('-ec', '--expressing-cell-count', default="num", type=str, help="column containing cluster id value")
    parser.add_argument('-tc', '--cluster-cell-count', default="anum", type=str, help="column containing cluster id value")

    parser.add_argument('-n', '--predictions', default=10, type=int, help="number of predictions per cluster shown")
    parser.add_argument('-f', '--mean-factor', default=1, type=float, help="number of predictions per cluster shown")
    parser.add_argument('-pf', '--pval-factor', default=1, type=float, help="number of predictions per cluster shown")

    parser.add_argument('-cdb', '--cellmarkerdb', default=False, action="store_true", help="generate seurat output at the end?")


    parser.add_argument('-s', '--seurat', default=False, action="store_true", help="generate seurat output at the end?")
    parser.add_argument('-sc', '--scanpy', default=False, action="store_true", help="generate scanpy output at the end?")
    
    parser.add_argument('-a', '--aorta3d', default=False, action="store_true", help="generate seurat output at the end?")
    parser.add_argument('-o', '--organs', default=[], type=str, nargs='+', help="generate seurat output at the end?")

    parser.add_argument('--output', type=argparse.FileType('w'), default=sys.stdout, help="write to output file")
    
    args = parser.parse_args()

    if args.seurat or args.aorta3d:
        print("Setting number of predictions to 1", file=sys.stderr)
        args.predictions = 1

    args.organs = [x.lower() for x in args.organs]


    print("Taking value gene from {}".format(args.gene), file=sys.stderr)
    print("Taking value cluster from {}".format(args.cluster), file=sys.stderr)
    print("Taking value logfc from {}".format(args.logfc), file=sys.stderr)
    print("Taking value pvaladj from {}".format(args.pvaladj), file=sys.stderr)
    print("Taking value expr-mean from {}".format(args.expr_mean), file=sys.stderr)
    print("Taking value expressing-cell-count from {}".format(args.expressing_cell_count), file=sys.stderr)
    print("Taking value cluster-cell-count from {}".format(args.cluster_cell_count), file=sys.stderr)

    gene2clusters = defaultdict(set)

    measuredGenes = set()

    with args.markers as fin:

        cluster2genes = defaultdict(lambda: dict())

        elem2idx = {}
        idx2elem = {}
        for idx, line in enumerate(fin):

            line = line.strip().split("\t")

            if idx == 0:

                for lidx, elem in enumerate(line):

                    elem2idx[elem] = lidx
                    idx2elem[lidx] = elem

                assert(args.gene in elem2idx)
                assert(args.cluster in elem2idx)
                assert(args.logfc in elem2idx)
                assert(args.pvaladj in elem2idx)
                assert(args.expr_mean in elem2idx)
                assert(args.expressing_cell_count in elem2idx)
                assert(args.cluster_cell_count in elem2idx)


                continue

            if line[0].startswith(idx2elem[0]):
                continue

            if len(line) < len(idx2elem):
                continue

            geneSym = line[elem2idx[ args.gene ]].upper()
            cluster = line[elem2idx[ args.cluster]]

            try:
                robustLogFC = 0
                if args.logfc != None:
                    robustLogFC = float(line[elem2idx[args.logfc]])
                robustPval = 0
                if args.pvaladj != None:
                    robustPval = float(line[elem2idx[args.pvaladj]]) * args.pval_factor
                expr_value = 0
                if args.expr_mean != None:
                    expr_value = float(line[elem2idx[args.expr_mean]]) * args.mean_factor
                expr_perc = 0
                if args.expressing_cell_count != None and args.cluster_cell_count != None:
                    expr_perc = float(line[elem2idx[args.expressing_cell_count]]) / float(line[elem2idx[args.cluster_cell_count]])

                if robustPval > 0.05:
                    continue

                measuredGenes.add(geneSym)

                # expr_value, logFC, pVal, _, epxr_perc
                cluster2genes[cluster][geneSym] = (expr_value, robustLogFC, robustPval, len(cluster2genes[cluster])+1, expr_perc)
            except:
                pass
        

    print("Got {} clusters.".format(len(cluster2genes)), file=sys.stderr)



    allFirstHits = {}

    if args.update_cellmarkerdb or not os.path.isfile("cellmarkerdb.tsv"):

        print("Updating CellMarker DB", file=sys.stderr)        

        url = "http://biocc.hrbmu.edu.cn/CellMarker/download/all_cell_markers.txt"

        try:
            with urllib.request.urlopen(url) as dl_file:
                with open("cellmarkerdb.tsv", 'wb') as outfile:
                    outfile.write(dl_file.read())
        except:
            print("Unable to download cellmarerdb", file=sys.stderr)

    if args.update_panglao or not os.path.isfile("panglao.tsv"):

        print("Did not find panglao file. Downloading it now", file=sys.stderr)

        try:

            url = "https://panglaodb.se/markers.html?cell_type=%27all_cells%27"
            user_agent = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:75.0) Gecko/20100101 Firefox/75.0'
            accept_header= "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8"
            request = urllib.request.Request(url,headers={'User-Agent': user_agent, 'accept': accept_header})
            response = urllib.request.urlopen(request)
            html = response.read()
            
            htmldoc = lxml.html.fromstring(html)
            allLinks = [x.attrib.get("href", None) for x in htmldoc.xpath(".//small/a")]

            allLinks = [x for x in allLinks if not x == None]

            for link in allLinks:
                if link.startswith("markers/"):
                    filelink = "https://panglaodb.se/{}".format(link)
                    print("Downloading from: ", filelink, file=sys.stderr)

                    if link.endswith(".gz"):
                        print("in compressed format", file=sys.stderr)

                        with urllib.request.urlopen(filelink) as dl_file:
                            with gzip.open(dl_file, 'rb') as fin:
                                with open("panglao.tsv", 'wb') as outfile:
                                    outfile.write(fin.read())

                    else:
                        print("in uncompressed format", file=sys.stderr)
                        with urllib.request.urlopen(filelink) as dl_file:
                            with open("panglao.tsv", 'wb') as out_file:
                                out_file.write(dl_file.read())

                break

        except:
            print("Unable to download panglao", file=sys.stderr)

    print("Starting analysis", file=sys.stderr)

    clusterid2genes = defaultdict(lambda: dict())
    gene2refcluster = defaultdict(set)
    gene2celltypes = defaultdict(set)
    nickname2gene = {}

    if not args.cellmarkerdb:

        with open("panglao.tsv") as fin:


            elem2idx = {}
            for idx, line in enumerate(fin):

                line = line.strip().split("\t")

                if idx == 0:

                    for lidx, elem in enumerate(line):
                        elem2idx[elem] = lidx

                    continue

                geneSym = line[elem2idx["official gene symbol"]]
                celltype = line[elem2idx["cell type"]]
                organ = line[elem2idx["organ"]]

                canonicalMarker = line[elem2idx["canonical marker"]] == "1"

                try:
                    meanSens = max(float(line[elem2idx["sensitivity_human"]]), float(line[elem2idx["sensitivity_mouse"]]))# / 2.0
                except:
                    
                    if line[elem2idx["sensitivity_human"]] != "NA":
                        meanSens = float(line[elem2idx["sensitivity_human"]])
                    elif line[elem2idx["sensitivity_mouse"]] != "NA":
                        meanSens = float(line[elem2idx["sensitivity_mouse"]])
                    else:
                        meanSens = 0

                try:
                    meanSpec = max(float(line[elem2idx["specificity_human"]]) + float(line[elem2idx["specificity_mouse"]]))# / 2.0
                except:
                    
                    if line[elem2idx["specificity_human"]] != "NA":
                        meanSpec = float(line[elem2idx["specificity_human"]])
                    elif line[elem2idx["specificity_mouse"]] != "NA":
                        meanSpec = float(line[elem2idx["specificity_mouse"]])
                    else:
                        meanSpec = 0

                measuredGenes

                nicknames = line[elem2idx["nicknames"]]
                if len(nicknames) > 0:

                    nicknames = nicknames.split("|")
                    replacedGeneSym = False

                    for nickn in nicknames:
                        nickname2gene[nickn] = geneSym

                        if not replacedGeneSym and not geneSym in measuredGenes:
                            if nickn in measuredGenes:
                                geneSym = nickn
                                replacedGeneSym = True

                if not geneSym in measuredGenes:
                    continue

                if canonicalMarker:
                    meanSens = 1.0

                gene2refcluster[geneSym].add((celltype, organ))

                clusterid2genes[(celltype, organ)][geneSym] = ( meanSens, meanSpec)
                gene2celltypes[geneSym].add(celltype)



    elif args.cellmarkerdb:
        print("Loading cellmarkerdb", file=sys.stderr)

        clusterid2genes = defaultdict(lambda: dict())

        tisCt2genes = defaultdict(set)
        gene2tisCt = defaultdict(set)

        with open("cellmarkerdb.tsv") as fin:

            elem2idx = {}
            for idx, line in enumerate(fin):

                line = line.strip().split("\t")

                if idx == 0:

                    for lidx, elem in enumerate(line):
                        elem2idx[elem] = lidx

                    continue

                #tissueType
                organ = line[elem2idx["tissueType"]]
                #cancerType
                isCancerType = not "NORMAL" in line[elem2idx["cellType"]].upper()
                #cellType
                #cellName
                celltype = line[elem2idx["cellName"]]

                if isCancerType:
                    continue

                if "et al" in celltype:
                    continue

                #cellMarker
                cellMarkers = line[elem2idx["cellMarker"]].replace("[", "").replace("]", "").split(", ")
                #geneSymbol
                geneSymbol = line[elem2idx["geneSymbol"]].replace("[", "").replace("]", "").split(", ")

                proteinName = line[elem2idx["proteinName"]].replace("[", "").replace("]", "").split(", ")

                tisCt = (celltype, organ)
                
                for gene in cellMarkers+geneSymbol+proteinName:
                    gene = gene.upper()

                    gene2tisCt[gene].add(tisCt)
                    tisCt2genes[tisCt].add(gene)

        with open("cellmarkerdb.tsv") as fin:

            clusterid2genes = defaultdict(lambda: dict())

            elem2idx = {}
            for idx, line in enumerate(fin):

                line = line.strip().split("\t")

                if idx == 0:

                    for lidx, elem in enumerate(line):
                        elem2idx[elem] = lidx

                    continue

                #tissueType
                organ = line[elem2idx["tissueType"]]
                #cancerType
                isCancerType = not "NORMAL" in line[elem2idx["cellType"]].upper()
                #cellType
                #cellName
                celltype = line[elem2idx["cellName"]]

                if isCancerType:
                    continue

                if "et al" in celltype:
                    continue

                #cellMarker
                cellMarkers = line[elem2idx["cellMarker"]].replace("[", "").replace("]", "").split(", ")
                #geneSymbol
                geneSymbol = line[elem2idx["geneSymbol"]].replace("[", "").replace("]", "").split(", ")

                tisCt = (celltype, organ)

                #if "MACROPHAGE" in celltype.upper():
                #    print(tisCt, tisCt2genes[tisCt])
               
                meanSpec = 1.0/len(tisCt2genes[tisCt])

                for gene in cellMarkers+geneSymbol+proteinName:
                    gene = gene.upper()

                    if not gene in measuredGenes:
                        continue

                    gene2refcluster[gene].add( tisCt )
                    clusterid2genes[tisCt][gene] = ( 1.0, 0.0)
                    gene2celltypes[gene].add(celltype)

 
    print("Loaded Databases", file=sys.stderr)
    print("known genes", len(gene2refcluster), file=sys.stderr)
    print("known (celltype, organ)", len(clusterid2genes), file=sys.stderr)



    for cluster in natsorted([x for x in cluster2genes]):
        clusterGenes = [x for x in cluster2genes[cluster]]

        clusterCounter = Counter()
        cluster2accGenes = dict()
        cluster2setAccGenes = dict()

        for refcluster in clusterid2genes:

            refclusterGenes = [x for x in clusterid2genes[refcluster]]

            totalScore = 0
            accGenes = 0
            setAccGenes = set()

            seenGenes = set()

            for gene in clusterGenes:

                geneFound = False
                if gene in refclusterGenes:
                    geneFound = True

                else:

                    if gene in nickname2gene:
                        gene = nickname2gene[gene]

                if gene in seenGenes:
                    continue

                seenGenes.add(gene)
                    
                if geneFound:
                    accGenes += 1
                    setAccGenes.add(gene)

                    # expr_value, logFC, pVal, rank, epxr_perc
                    geneInfo = cluster2genes[cluster][gene]
                    clusterInfo = clusterid2genes[refcluster][gene]
                    #sens, spec

                    importanceInRefCluster = 1.0 / (1+math.log2(len(gene2refcluster[gene])))

                    importanceForCellType = 1.0

                    if len(gene2celltypes[gene]) == 1:
                        importanceForCellType = 1.5

                    meanSens = clusterInfo[0]
                    meanSpec = clusterInfo[1]

                    pvalImp = 1000

                    if geneInfo[2] == 0.0:
                        pvalImp = 1000
                    else:
                        pvalImp = -math.log(geneInfo[2], 2)


                    avgExpr = geneInfo[0]
                    logFC = geneInfo[1]
                    geneRank = geneInfo[3]
                    prevInCluster = geneInfo[4]

                    fGeneRank = 1 / math.log(geneRank+1, 2)

                    if logFC > 0 and prevInCluster > 0.5:
                        #was  meanSens * logFC * prevInCluster
                        #geneScore = meanSens * logFC  * (1.0-meanSpec) * importanceInRefCluster # v1.0
                        #geneScore = meanSens * logFC  * (1.0-meanSpec) * importanceInRefCluster * importanceForCellType # version 1.1
                        geneScore = meanSens * avgExpr * prevInCluster * (1.0-meanSpec) * importanceInRefCluster  # version 1.2
                        #print(gene, geneScore)
                        totalScore += geneScore

                    #else:
                    #    print(gene, logFC)
                        

                cluster2accGenes[refcluster] = accGenes
                cluster2setAccGenes[refcluster] = setAccGenes

                accGenesUniqueForCelltype = 0
                genesUniqueForCelltype = 0

                ctUniqueGenes = set()
                for g in clusterid2genes[refcluster]:

                    if len(gene2celltypes[g]) == 1:
                        genesUniqueForCelltype += 1

                        if g in cluster2setAccGenes[refcluster]:
                            accGenesUniqueForCelltype += 1
                            ctUniqueGenes.add(g)

            if genesUniqueForCelltype == 0:
                accGenesUniqueForCelltype = 1
                genesUniqueForCelltype = 10
            elif accGenesUniqueForCelltype == 0:
                accGenesUniqueForCelltype = 0.1

            if accGenes == 0:
                continue
            #clusterCounter[refcluster] = totalScore * (accGenes/len(clusterid2genes[refcluster])) # len(clusterid2genes[refcluster]) v1.0
            clusterCounter[refcluster] = totalScore * (accGenesUniqueForCelltype/genesUniqueForCelltype) * (accGenes/len(clusterid2genes[refcluster])) # len(clusterid2genes[refcluster]) v1.1


        accOutput = 0
        for idx, x in enumerate(clusterCounter.most_common()):

            if accOutput < args.predictions:
                
                if args.organs != None and len(args.organs) > 0:

                    thisOrgan = x[0][1].lower()
                    if not thisOrgan in args.organs:
                        continue
                
                accGenesUniqueForCelltype = 0
                genesUniqueForCelltype = 0

                ctUniqueGenes = set()
                for g in clusterid2genes[x[0]]:

                    if len(gene2celltypes[g]) == 1:
                        genesUniqueForCelltype += 1

                        if g in cluster2setAccGenes[x[0]]:
                            accGenesUniqueForCelltype += 1
                            ctUniqueGenes.add(g)

                print(cluster, ";".join(x[0]), x[1], cluster2accGenes[x[0]], len(clusterid2genes[x[0]]),
                accGenesUniqueForCelltype, genesUniqueForCelltype, ctUniqueGenes, cluster2setAccGenes[x[0]], sep="\t", file=args.output)
                if accOutput == 0:
                    allFirstHits[cluster] = ";".join(x[0])
                accOutput += 1


    if args.seurat:

        outstr = "new.cluster.ids <- c({})".format(
            ",".join(['"{}"'.format(allFirstHits[x]) for x in allFirstHits])
        )

        print(outstr)
        print("orignames = Idents(seurat_obj)")
        print("names(new.cluster.ids) <- levels(orignames)")
        print("levels(orignames) = new.cluster.ids")
        print("seurat_obj$cellnames = orignames")
        
    elif args.scanpy:
        
        outstr = "group2cellname <- dict({})".format(
            ",".join(['"{}" = "{}"'.format(x, allFirstHits[x]) for x in allFirstHits])
        )

        print(outstr)
        
        scanpycode = """
{o2ndict}  
  
adata.obs['new_clusters'] = (
    adata.obs['leiden_0.6']
    .map(group2cellname)
    .astype('category')
)
        """.format(o2ndict=outstr)

        print(scanpycode)
        
