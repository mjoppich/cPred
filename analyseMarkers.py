from collections import defaultdict, Counter
import math
import sys, os
import argparse
import gzip
from natsort import natsorted

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Name Seurat clusters')
    parser.add_argument('-i', '--markers', type=argparse.FileType('r'),help='output from Seurat FindMarkers or FindAllMarkers')
    parser.add_argument('-u', '--update-panglao', action='store_true', default=False, help='update panglao db file')

    parser.add_argument('-g', '--gene', default="gene", type=str, help="column containing cluster id value")
    parser.add_argument('-c', '--cluster', default="clusterID", type=str, help="column containing cluster id value")

    parser.add_argument('-l', '--logfc', default="avg_logFC", type=str, help="column containing cluster id value")
    parser.add_argument('-p', '--pvaladj', default="p_val_adj", type=str, help="column containing cluster id value")
    parser.add_argument('-e', '--expr-mean', default="mean", type=str, help="column containing cluster id value")

    parser.add_argument('-ec', '--expressing-cell-count', default="num", type=str, help="column containing cluster id value")
    parser.add_argument('-tc', '--cluster-cell-count', default="anum", type=str, help="column containing cluster id value")

    parser.add_argument('-n', '--predictions', default=10, type=int, help="number of predictions per cluster shown")

    parser.add_argument('-s', '--seurat', default=False, action="store_true", help="generate seurat output at the end?")
    
    args = parser.parse_args()

    if args.seurat:
        print("Setting number of predictions to 1", file=sys.stderr)
        args.predictions = 1



    gene2clusters = defaultdict(set)

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

                continue

            if line[0].startswith(idx2elem[0]):
                continue

            if len(line) < len(idx2elem):
                continue

            geneSym = line[elem2idx[ args.gene ]]
            cluster = line[elem2idx[ args.cluster]]

            try:
                robustLogFC = 0
                if args.logfc != None:
                    robustLogFC = float(line[elem2idx[args.logfc]])
                robustPval = 0
                if args.logfc != None:
                    robustPval = float(line[elem2idx[args.pvaladj]])
                expr_value = 0
                if args.logfc != None:
                    expr_value = float(line[elem2idx[args.expr_mean]])
                expr_perc = 0
                if args.logfc != None:
                    expr_perc = float(line[elem2idx[args.expressing_cell_count]]) / float(line[elem2idx[args.cluster_cell_count]])

                if robustPval > 0.05:
                    continue

                cluster2genes[cluster][geneSym] = (expr_value, robustLogFC, robustPval, len(cluster2genes[cluster])+1, expr_perc)
            except:
                pass
        

    gene2refcluster = defaultdict(set)
    nickname2gene = {}

    allFirstHits = []

    if args.update_panglao or not os.path.isfile("panglao.tsv"):

        print("Did not find panglao file. Downloading it now", file=sys.stderr)

        import urllib.request
        import lxml.html

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

    print("Starting analysis", file=sys.stderr)

    with open("panglao.tsv") as fin:

        clusterid2genes = defaultdict(lambda: dict())

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


            gene2refcluster[geneSym].add((celltype, organ))

            clusterid2genes[(celltype, organ)][geneSym] = ( meanSens, meanSpec)

            nicknames = line[elem2idx["nicknames"]]
            if len(nicknames) > 0:

                nicknames = nicknames.split("|")

                for nickn in nicknames:
                    nickname2gene[nickn] = geneSym



    for cluster in natsorted([x for x in cluster2genes]):

        clusterGenes = [x for x in cluster2genes[cluster]]

        clusterCounter = Counter()
        cluster2accGenes = dict()

        for refcluster in clusterid2genes:

            refclusterGenes = [x for x in clusterid2genes[refcluster]]

            totalScore = 0
            accGenes = 0

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
                    geneInfo = cluster2genes[cluster][gene]
                    clusterInfo = clusterid2genes[refcluster][gene]
                    #sens, spec

                    importanceInRefCluster = 1.0 / (1+math.log2(len(gene2refcluster[gene])))

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

                    if logFC > 0:
                        #was  meanSens * logFC * prevInCluster
                        totalScore += meanSens * logFC  * (1.0-meanSpec) * importanceInRefCluster #(1-meanSpec) *fGeneRank *pvalImp
                        
            clusterCounter[refcluster] = totalScore * (accGenes/len(clusterid2genes[refcluster])) # len(clusterid2genes[refcluster])
            cluster2accGenes[refcluster] = accGenes

        for idx, x in enumerate(clusterCounter.most_common(args.predictions)):
            print(cluster, ";".join(x[0]), x[1], cluster2accGenes[x[0]], len(clusterid2genes[x[0]]), sep="\t")

            if idx == 0:
                allFirstHits.append(";".join(x[0]))

    if args.seurat:

        outstr = "new.cluster.ids <- c({})".format(
            ",".join(['"{}"'.format(x) for x in allFirstHits])
        )

        print(outstr)
        print("orignames = Idents(seurat_obj)")
        print("names(new.cluster.ids) <- levels(orignames)")
        print("levels(orignames) = new.cluster.ids")
        print("seurat_obj$cellnames = orignames")
