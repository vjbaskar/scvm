from .globimport import *

class Gsea:
    """
    Run GSEA.
    Uses Enrichr databases
    
    adata: anndata
    groupby: obs column name based on which you have run sc.tl.rank_gene_groups()
    organism: name of organism
    gene_sets: gene_sets to use. Run Gsea.getdb() to get gene lists available.
    de_pval_cutoff: pval_cutoff to use for rank_gene_groups results to get genes
    de_log2fc_min: min log2fc to be for rank_gene_groups results to get genes
    top_genes: [int] top genes to be used.
    outdir: Store results in a dir.
    
    https://gseapy.readthedocs.io/en/latest/introduction.html
    
    In adata the results are stored in adata.uns['gsea']
    
    Example:
    gsea = cf.Gsea(adata)
    gsea.run_gsea()
    adata = gsea.get_results()
    gsea.get_results(group='14', pcutoff=0.05)
    
    """
    def __init__(self, adata, groupby = 'leiden', organism = 'mouse'):
        self.adata = adata
        self.groupby = groupby
        self.adata.uns['gsea'] = dict()
        self.organism = organism
    
    def _gsea(self, i, gene_sets, organism, **kwargs):
        """
        Internal function to run GSEA
        """
        import gseapy as gp

        adata = self.adata
        top_genes = self.top_genes
        print(f"Running Enrichr for group:{i} through {gene_sets}")
        
        df = sc.get.rank_genes_groups_df(adata, group=i, pval_cutoff=self.de_pval_cutoff, log2fc_min=self.de_log2fc_min)
        
        if df.shape[0] < top_genes:
            top_genes = df.shape[0]
        print(f"total genes that passed cutoff = {df.shape[0]} | genes taken is {top_genes}")
        
        genelist = df.names[0:100].str.upper().values.tolist()
        enr = gp.enrichr(gene_list=genelist, # or "./tests/data/gene_list.txt",
                     gene_sets=['PanglaoDB_Augmented_2021'],
                     organism=organism, 
                             **kwargs
                    )
        return enr.results
    
    def getdb(self):
        """
        Get available databases
        """
        x = gp.get_library_name(organism=self.organism)
        from pprint import pprint as pp
        pp(x)
    
    def run_gsea(self, gene_sets=['PanglaoDB_Augmented_2021'], group = None, de_pval_cutoff = 0.05, de_log2fc_min = 1, top_genes = 100, outdir = None, **kwargs):
        
        """
        Run GSEA
        """
        
        self.de_pval_cutoff = de_pval_cutoff
        self.de_log2fc_min = de_log2fc_min
        self.top_genes = top_genes
        
        adata = self.adata
        if group is None:
            for i in self.adata.obs[self.groupby].unique():
                self.adata.uns['gsea'][i] = self._gsea(i, gene_sets, organism = self.organism, outdir = None)
        else:
            i = group
            self.adata.uns['gsea'][i] = self._gsea(i, gene_sets, organism = self.organism, outdir = None)
        
    def get_results(self, adata = None, group = None, pcutoff = 0.05):
        if adata is None:
            temp = self.adata
        else:
            temp = adata
        if group is None:
            return self.adata
        else:
            x = temp.uns['gsea'][group]
            x = x[x['Adjusted P-value'] <= pcutoff]
            return x



