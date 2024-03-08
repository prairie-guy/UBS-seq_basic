<<CBD: How are each of these used? >>

CUSTOMIZED_GENES = [] <<CBD: Are these genes to be used in experiment; Need reference fasta? >>
Filter out tRNA or rRNA in mapping genes  hisat2_3n_mapping _genes
Not need for DNA

LIBRARY_STRATEGY = "INLINE" <<CBD: What is INLINE? Other strategies: TAKARAV3, SWIFT, STRANDED. Why and how are each used? >>
INLINE is the default. STRANDED is a reverse complement (No UMI) R

STRANDNESS = "F" <<CBD: This is the default. Where is this relevent? Cutadapt, alignment? >>
INLINE based upon F

REFTYPES = ["genes", "genome", "contamination"] <<CBD: Is this used for reference fasta? What is the strategy here? I did not see any filtering of rRNA, tRNA or contamination? Are these run for each sample? >>
Mapping: 1 mapping genes, 2 mapping genoome 3 mapping contaminaiton With each feeding the next.

REFTYPES_CALL = ["genes", "genome"] <<CBD: Is genes for mRNA and genome for DNA? How are these set in config.yaml? >>

pairend_run_ids = [] <<CBD: What are these used for? >>
Not major

for s, v in config["samples"].items(): <<CBD: Are groups considered to be replicates? Are runs elements of a group? >>
Samples and groups are major concepts

<<CBD: How far through the pipeline must paired ends be condidered seperately> qc, cutadapt, hisat-3n_mapping, sorting, joining >>
Come together at calling

<<CBD: Prepared for customized. For default which version of index used? --repeat-index required 256 Gb RAM. What does this do? Can I use yours. Does any version of GRCh38 work? Using v.110 from Ensemble. >>

As long as major are same, ok to use different annotations and indexes.

<<CBD: What is directional mapping? --directional-mapping (F) vs. --directional-mapping-reverse (R) vs. nothing>>
DNA are forward. Only relevant RNA

<<CBD: What does --norc (No reverse complement) >>
For gene mapping, not required. Only for genome. Mostly optimization related to avoid 2x mapping

<<CBD: How to use --no-spliced-alignment >>
Not needed for genes, but needed for genome for RNA data. (For DNA not needed at all)

<<CBD: In joining SE and PE mappings, is it typical to have SE joined with PE? What is being merged? It looks like runs (same as replicates) are kept seperate >>
In some cases it might happen. Rare. Corner case. Ignore

<<CBD: In combining runs, is this the same as merging all groups together. For example, all of treated and all of control. For example, does t1,t2 -> t and c1 -> c >>
Yes, but still have tables for user for downstream analysis.
Runs are a differnt concepct than groups. Run = Seq Run vs. Biological replicate.

<<CBD: In combining runs, do we lose the ability for to look at replicates against one another>>

<<CBD: Need parse_cutadapt_report.py>>
Done

<<CBD: What is being done here? What is the output>>

<<CBD: stat_mapping looks to be creating a tsv file for REFTYPES = ["genes", "genome", "contamination"]. Finds number  of high quality reads>>


<<CBD: Why not filter out these reads? >>

<<CBD: Does dedup only work on combined groups. Are single samples not deduped?>>
Combine and then dedup seq runs

<<CBD: Call mutations for unique reads from deduped combined reads.  What is importance of unique vs. multi aligned >>


<<CBD: Call mutations for multi reads from deduped combined reads>>

<<CBD: Why filter -e "rlen<100000" >>

<<CBD: Filter (uncalled) combined deduped >>

<<CBD: Why only < 4 mutations? >>
Because its on a read basis

<<CBD: Technically output is not really converted>>

<<CBD: Call unique previously filtered, combined >>

<<CBD: Call multiple aligned,  previously filtered, combined >>

<<CBD: How to get CG context. CGH etc >>
Bismark does this, but would have to be done in parallel
htslib is C library for samtools
Python script

<<CBD: Is this being called on all called samples? >>

<<CBD: Appears to be nuc/(nuc + nc)

<<CBD: Need select_called_sites_v3.py >>

<<CBD: What is logic? >>

<<CBD: Need combined_selected_sites.py >>

<<CBD: What is logic? >>

Background
Built to determine the background
