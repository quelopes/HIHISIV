sql_max_abs_log_fc = """
SELECT      round(max(abs(log_fc)),1)-0.1
FROM        deg_analysis 
WHERE       hgnc_symbol!='';
"""

sql_list_genes_by_symbol = """
SELECT      DISTINCT(deg.hgnc_symbol)
FROM        public.deg_analysis deg 
LEFT JOIN   experiment expe ON expe.experiment_id = deg.experiment_id 
LEFT JOIN   experimental_factor_test test ON test.experiment_id = expe.experiment_id 
LEFT JOIN   experimental_factor_reference refe ON refe.experiment_id = expe.experiment_id 
WHERE       abs(deg.log_fc) >= %(logFC)s AND deg.adj_pvalue<=%(adjPvalue)s AND 
            deg.hgnc_symbol != '' order by deg.hgnc_symbol
"""

sql_list_genes_by_entrez_id = """
SELECT      DISTINCT(deg.entrez_id)
FROM        public.deg_analysis deg 
LEFT JOIN   experiment expe ON expe.experiment_id = deg.experiment_id 
LEFT JOIN   experimental_factor_test test ON test.experiment_id = expe.experiment_id 
LEFT JOIN   experimental_factor_reference refe ON refe.experiment_id = expe.experiment_id 
WHERE       abs(deg.log_fc) >= %(logFC)s AND deg.adj_pvalue<=%(adjPvalue)s AND 
            deg.entrez_id != '' order by deg.entrez_id
"""

sql_gene_expression_experiment_by_symbol = """
SELECT      DISTINCT expe.experiment_id as experiment_id, 
            expe.sample_description as experiment_sample_description,
            deg.hgnc_symbol as symbol,
            to_char(log_fc, '999D99') as log_fc, 
            to_char(adj_pvalue, '9.99EEEE') as adj_pvalue, 
            n_samples,
            reference_factor_name as reference_name, 
            reference_factor_ontology_id as reference_ontology, 
            reference_sp as reference_species, 
            reference_virus,
            test_factor_name as test_name, 
            test_factor_ontology_id as test_ontology, 
            test_sp as test_species,
            test_virus
FROM        public.deg_analysis deg 
LEFT JOIN   experiment expe ON expe.experiment_id = deg.experiment_id
LEFT JOIN   experimental_factor_test test ON test.experiment_id = expe.experiment_id 
LEFT JOIN   experimental_factor_reference refe ON refe.experiment_id = expe.experiment_id 
WHERE       deg.hgnc_symbol = %(gene_sel_id)s AND abs(deg.log_fc) >= %(logFC)s and deg.adj_pvalue<=%(adjPvalue)s;
"""

sql_gene_expression_experiment_by_entrez_id = """
SELECT      DISTINCT expe.experiment_id as experiment_id, 
            expe.sample_description as experiment_sample_description,
            deg.hgnc_symbol as symbol,
            to_char(log_fc, '999D99') as log_fc, 
            to_char(adj_pvalue, '9.99EEEE') as adj_pvalue, 
            n_samples,
            reference_factor_name as reference_name, 
            reference_factor_ontology_id as reference_ontology, 
            reference_sp as reference_species, 
            reference_tissue_ontology as reference_tissue, 
            reference_virus,
            test_factor_name as test_name, 
            test_factor_ontology_id as test_ontology, 
            test_sp as test_species,
            test_tissue_ontology as test_tissue, 
            test_virus 
FROM        public.deg_analysis deg 
LEFT JOIN   experiment expe ON expe.experiment_id = deg.experiment_id
LEFT JOIN   experimental_factor_test test ON test.experiment_id = expe.experiment_id 
LEFT JOIN   experimental_factor_reference refe ON refe.experiment_id = expe.experiment_id 
WHERE       deg.entrez_id = %(gene_sel_id)s AND abs(deg.log_fc) >= %(logFC)s and deg.adj_pvalue<=%(adjPvalue)s;
"""

sql_list_transcripts =  """
SELECT      DISTINCT(transcript_id)
FROM        public.deg_analysis deg
LEFT JOIN   experiment expe ON expe.experiment_id = deg.experiment_id
LEFT JOIN   experimental_factor_test test ON test.experiment_id = expe.experiment_id
LEFT JOIN   experimental_factor_reference refe ON refe.experiment_id = expe.experiment_id
WHERE       abs(deg.log_fc) >= %(logFC)s AND deg.adj_pvalue<=%(adjPvalue)s AND
            transcript_id != '' order by transcript_id
"""

sql_transcript_experiment = """
SELECT      DISTINCT expe.experiment_id as experiment_id,
            expe.sample_description as experiment_sample_description,
            deg.hgnc_symbol as symbol,
            to_char(log_fc, '999D99') as log_fc, 
            to_char(adj_pvalue, '9.99EEEE') as adj_pvalue, 
            n_samples,
            reference_factor_name as reference_name,
            reference_factor_ontology_id as reference_ontology,
            reference_sp as reference_species,
            reference_tissue_ontology as reference_tissue,
            reference_virus,
            test_factor_name as test_name,
            test_factor_ontology_id as test_ontology,
            test_sp as test_species,
            test_tissue_ontology as test_tissue,
            test_virus
FROM        public.deg_analysis deg
LEFT JOIN   experiment expe ON expe.experiment_id = deg.experiment_id
LEFT JOIN   experimental_factor_test test ON test.experiment_id = expe.experiment_id
LEFT JOIN   experimental_factor_reference refe ON refe.experiment_id = expe.experiment_id
WHERE       deg.transcript_id = %(transcript_sel_id)s AND abs(deg.log_fc) >= %(logFC)s AND
            deg.adj_pvalue<=%(adjPvalue)s
"""

sql_list_go_terms =  """
SELECT      DISTINCT(go_term)
FROM        go_analysis 
ORDER BY    go_term
"""

sql_get_go_id = """
SELECT      go_id
FROM        go_analysis
WHERE       go_term=%(go_term)s;
"""

sql_experiment_go_term =  """
SELECT      DISTINCT goa.experiment_id, 
            goa.go_id, 
            goa.go_term, 
            goa.p_value,
            goa.no_go_size,
            expe.sample_description,	
            test.test_factor_name, 
            refe.reference_factor_name, 
            goa.entrez 
FROM        go_analysis goa 
LEFT JOIN   deg_analysis deg ON deg.experiment_id = goa.experiment_id 
LEFT JOIN   experiment expe ON expe.experiment_id = goa.experiment_id 
LEFT JOIN   experimental_factor_test test ON test.experiment_id = goa.experiment_id 
LEFT JOIN   experimental_factor_reference refe ON refe.experiment_id = goa.experiment_id 
WHERE       go_term=%(go_term)s
"""

sql_list_ontology_terms = """
SELECT      DISTINCT(term) 
FROM        ontology
WHERE       term!=''
ORDER BY    term;
""" 

sql_experiment_ontology_terms = """
SELECT          DISTINCT db_accession, 
		onto.term, 		
		onto.ontology_name,
		onto.ontology_id,
		onto.tag as ontology_tag,
        experiment.sample_description as experiment_sample_description,
		plat.geo_platform_id as platform_id, 
		proj.title,
		proj.summary
FROM            ontology onto NATURAL JOIN experiment
NATURAL JOIN    project as proj
NATURAL JOIN    platform_project plat
WHERE           onto.term 
LIKE %(term_sel)s
"""

sql_single_gene_co_expression_graph =  """
WITH experiment_gene AS (
    SELECT DISTINCT experiment_id, hgnc_symbol 
    FROM            deg_analysis 
    WHERE           hgnc_symbol!='' AND adj_pvalue<=%(adjPvalue)s and abs(log_fc)>=%(logFC)s
    )
SELECT      hgnc_symbol, count(*) as value, string_agg(experiment_id, ', ') as experiment_list 
FROM        experiment_gene  
WHERE       hgnc_symbol!=%(gene_symbol)s AND experiment_id in (
                SELECT  experiment_id
                FROM    experiment_gene 
                WHERE   hgnc_symbol=%(gene_symbol)s
                )
GROUP BY    hgnc_symbol order by value desc limit 50
"""

sql_co_expressed_genes = """
WITH experiment_gene AS (
    SELECT DISTINCT experiment_id, hgnc_symbol 
    FROM            deg_analysis 
    WHERE           hgnc_symbol!='' AND adj_pvalue<=%(adjPvalue)s and abs(log_fc)>=%(logFC)s
    )
SELECT      hgnc_symbol, count(*) as value, string_agg(experiment_id, ', ') as experiment_list
FROM        experiment_gene  
WHERE       hgnc_symbol!=%(gene_symbol)s AND experiment_id in (
                SELECT  experiment_id
                FROM    experiment_gene 
                WHERE   hgnc_symbol=%(gene_symbol)s
                )
GROUP BY    hgnc_symbol
HAVING      count(*)>=%(co_threshold)s 
ORDER BY    value DESC
"""
