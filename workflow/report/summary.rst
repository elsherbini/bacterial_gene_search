Bacterial Gene Search: Summary Report
===============================

This report summarizes the results of searching for Pfam domains 
across the provided bacterial genomes.

Searched Pfams
-------------

The following Pfam domains were searched:

{% for row in snakemake.params.pfams %}
* **{{ row.pfam_accession }}**: {{ row.description }}
{% endfor %}

Searched Genomes
---------------

A total of {{ snakemake.params.sample_count }} genomes were searched:

{% for sample in snakemake.params.samples %}
* {{ sample }}
{% endfor %}

Results Summary
-------------

Here is a summary of the genes found in each genome:

{% for gene in snakemake.params.genes %}
### {{ gene.pfam_accession }}: {{ gene.description }}

Found in {{ gene.found_in|length }} out of {{ snakemake.params.sample_count }} genomes.

{% if gene.found_in|length > 0 %}
**Present in:**
{% for sample in gene.found_in %}
* {{ sample }}
{% endfor %}
{% endif %}

{% if gene.tree_available %}
[View phylogenetic tree]({{ gene.pfam_accession }}.aa.tree)
{% endif %}

{% endfor %} 