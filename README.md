# scrnaseq_celltype_prediction

Predict the cell-type of a cluster from Seurat FindMarker data

This method has not been tested on large scale data, but might actually work.

If you find this worked for you, please drop a line :)

## Usage

Download the analyseMarkers.py script and run it on your output files.

If you have a list with marker genes for each cluster (`list(clusterID0=..., clusterID1=...)`) then the best workflow to get expression values and a compatible data frame is: from within R/Seurat run the following script to generate a dataframe containing both marker gene values and expression data (code not yet tested ...):

    getExprData = function(markerObj, markerCells)
    {
    expTable = GetAssayData(object = subset(x=markerObj, cells=markerCells), slot = "data")
    
    outvalues1 = t(apply(expTable, 1, function(x) {
        a=x[x > 0];
        #a=x;
        out = {}
        
        out["anum"] = length(x)
        out["num"] = length(a)
        
        f = fivenum(a)
        out["min"] = f[1]
        out["lower_hinge"] = f[2]
        out["median"] = f[3]
        out["upper_hinge"] = f[4]
        out["max"] = f[5]
        out["mean"] = mean(a)
        
        out
    }))
    outvalues1 = cbind(rownames(outvalues1), outvalues1)
    cnames = colnames(outvalues1)
    cnames[1] = "gene"
    colnames(outvalues1) = cnames
    
    return(outvalues1)
    }

    getDEXpressionDF = function ( scdata, markers, assay="SCT" )
    {
    
    outDF = NULL
    DefaultAssay(object=scdata) = assay  
    clusterIDs = as.character(sort(unique(Idents(scdata))))
    
    scCells = Idents(scdata)
    scCells = names(scCells)
    scCells = unlist(as.character(scCells))
    
    for (clusterID in clusterIDs){
        
        print(clusterID)
        
        cellIdents = Idents(scdata)
        cellIdents.c = names(cellIdents[cellIdents == clusterID])
        cellIdents.c = unlist(lapply(cellIdents.c, as.character))    
        
        expvals = getExprData(scdata, cellIdents.c)

        modmarkers = markers[[clusterID]]
        modmarkers$gene = rownames(modmarkers)
        
        markerdf = as.data.frame(modmarkers)
        
        if ((nrow(markerdf) > 0) && (nrow(expvals) > 0))
        {
        expvals = merge(markerdf, expvals, all.x=T, by.x="gene", by.y = "gene")  
        }
        
        expvals = as.data.frame(cbind(clusterID, expvals))
        
        if (!is.data.frame(outDF) || nrow(outDF)==0)
        {
        outDF = expvals
        } else {
        outDF = as.data.frame(rbind(outDF, expvals))
        }
        
    }
    
    return(outDF)
    
    }

    makeDEResults = function(inobj, assay="SCT", test="wilcox")
    {
    clusterIDs = as.character(sort(unique(Idents(inobj))))

    retList = list()

    for (clusterID in clusterIDs)
    {


        cellIdents = Idents(inobj)
        cellIdents.c = names(cellIdents[cellIdents == clusterID])
        cellIdents.c = unlist(lapply(cellIdents.c, as.character))

        print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))

        deMarkers = FindMarkers(inobj, assay=assay, ident.1 = cellIdents.c, test.use=test)


        retList[[clusterID]] = deMarkers

    }

    return(retList)

    }


    deRes = makeDEResults(seurat_obj, assay="RNA", test="MAST")
    exprdf = getDEXpressionDF(seurat_obj, deRes, assay="RNA")
    write.table(exprdf, "expr_test.tsv", sep="\t", row.names=F, quote = F)



Then you can simply call

    python3 analyseMarkers.py --markers example/marker_genes_single_human_all.tsv --predictions 1
    Did not find panglao file. Downloading it now
    Downloading from:  https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz
    in compressed format
    Starting analysis
    0       Fibroblasts;Connective tissue   3.7979322878001667      62      179
    1       Smooth muscle cells;Smooth muscle       3.209013271632466       38      82
    2       Pancreatic stellate cells;Pancreas      2.2069340605035954      20      29
    3       Macrophages;Immune system       5.050242002842474       63      153
    4       Smooth muscle cells;Smooth muscle       0.9633089909462268      23      82
    5       Endothelial cells;Vasculature   7.505313125200414       84      195
    6       Fibroblasts;Connective tissue   2.8734895530493767      55      179
    7       T cells;Immune system   1.4604815939740858      35      107
    8       Fibroblasts;Connective tissue   0.7440482333103114      34      179
    9       Fibroblasts;Connective tissue   3.883677668408906       63      179
    10      Pericytes;Vasculature   0.803543270719957       22      64
    11      Smooth muscle cells;Smooth muscle       1.985453728430751       37      82
    12      B cells naive;Immune system     0.9207522036965257      19      69
    13      Plasma cells;Immune system      13.783575029154143      39      86
    14      Dendritic cells;Immune system   3.998198909943591       49      133
    15      Kupffer cells;Liver     1.9725414603905738      16      49
    16      Pancreatic stellate cells;Pancreas      3.883768523344435       18      29
    17      Endothelial cells;Vasculature   5.854839914806603       68      195
    18      Endothelial cells;Vasculature   6.300630147493595       78      195
    19      Schwann cells;Brain     2.201608826888472       21      48
    20      Fibroblasts;Connective tissue   0.5834923037164054      37      179
    21      NK cells;Immune system  3.665464137900602       40      98
    23      Mast cells;Immune system        1.9691749668080014      21      162

and you will receive 1 prediction per cluster. In a real life scenario you might want to check more than just one.

The output format is as follow:
    cluster -> cell_type -> score -> accepted_marker_genes -> marker_genes_of_celltype


### Renaming clusters in Seurat

With the `--seurat` flag, this tool will generate string that can easily be pasted into your R session:

    new.cluster.ids <- c("Fibroblasts;Connective tissue", "Smooth muscle cells;Smooth muscle", "Fibroblasts;Connective tissue", "Macrophages;Immune system", "Smooth muscle cells;Smooth muscle", "Endothelial cells;Vasculature", "Fibroblasts;Connective tissue", "T memory cells;Immune system", "Fibroblasts;Connective tissue", "Fibroblasts;Connective tissue", "Pericytes;Vasculature", "Smooth muscle cells;Smooth muscle", "B cells;Immune system", "Plasma cells;Immune system", "Macrophages;Immune system", "Macrophages;Immune system", "Fibroblasts;Connective tissue", "Endothelial cells;Vasculature", "Endothelial cells;Vasculature", "Schwann cells;Brain", "Fibroblasts;Connective tissue", "Gamma delta T cells;Immune system", "Mesothelial cells;Epithelium", "Mast cells;Immune system")

    orignames = Idents(seurat_obj)
    names(new.cluster.ids) <- levels(orignames)
    levels(orignames) = new.cluster.ids

    seurat_obj$cellnames = orignames

You can visualize the assigned cell types in a UMAP plot with

    DimPlot(obj.integrated, group.by="cellnames", reduction = "umap", label=T)

## Method

This prediction tools makes use of the marker genes provided by PanglaoDB [1]. Together with the reported sensitivity and specificity reported by them as well, the provided marker genes per cell-type and tissue are important. The script will download this marker table automatically.

In their original publication [1] Franzén et al. propose a weight sum for determining the cell type in scRNAseq experiments, however do not provide an implementation that could be executed on the data directly.

Here a slightly modified version of the score is implemented.

The prediction of the cell type for a specfic cluster is achieved using a weighted sum.
For each cell type $j$, and for each expressed gene *g* in a cluster *k**, which also happens to be a marker gene (it will be named *accepted* gene), a gene score 

![formula](https://latex.codecogs.com/png.latex?%5Ctext%7BGS%7D_%7Bj%2Ck%2Cg%7D%20%3D%20%5Ctext%7BmeanSens%7D_%7Bj%2Cg%7D%20%5Ccdot%20%5Ctext%7BavgExpr%7D_%7Bg%2Ck%7D%20%5Ccdot%20%5Ctext%7BprevInCluster%7D_%7Bk%2Cg%7D%5Ccdot%20%5Ctext%7BimpRC%7D_%7Bg%7D%20%5Ccdot%20%281-%5Ctext%7BmeanSpec%7D_%7Bj%2Cg%7D%29)


is calculated.
Here *meanSens* refers to the sensitivity with which the gene is expressed in the cell type, likewise *meanSpec* refers to the associated specificity.
Contributing to this score is also the prevalence of the gene in the cluster, *prevInCluster* - a measure that is deducted from the cluster annotation, namely the number of cells which express the gene and the total cells in the cluster.
The average expression *avgExpr* relates to the mean expression of the specific gene in the cluster.

This gene score is summed up for each accepted gene of a cluster such that the totalScore 

![formula](https://latex.codecogs.com/png.latex?%5Ctext%7BtotalScore%7D_%7Bj%2Ck%7D%20%3D%20%5Csum_%7Bg%20%5Cin%20MG_k%7D%20%5Ctext%7BGS%7D_%7Bj%2Ck%2Cg%7D)

for all significant marker genes *MG_k* of cluster *k*.

The final cluster score is then


![formula](https://latex.codecogs.com/png.latex?%5Ctext%7BclusterScore%7D_%7Bj%2Ck%7D%20%3D%20%5Ctext%7BtotalScore%7D_%7Bj%2Ck%7D%20%5Ccdot%20%5Cfrac%7B%5Ctext%7BaccUnique%7D_%7Bj%2Ck%7D%7D%7B%5Ctext%7BallUnique%7D_%7Bj%7D%7D%20%5Ccdot%20%5Cfrac%7B%5Ctext%7BaccGenes%7D_%7Bj%2Ck%7D%7D%7B%5Ctext%7BallGenes%7D_%7Bj%7D%7D)

where accUniquerefers to the in cluster *k* accepted unique genes for cell type *j*, and allUnique to all unique genes of that specific cell type.
Likewise accGenes is the number of accepted genes in cluster *k* for cell type *j*, and allGenes the number of all marker genes for cell type $j$.


[1] O. Franzén, L.-M. Gan, and J. L. M. Björkegren, “PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data,” Database, vol. 2019, Jan. 2019.