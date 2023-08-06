DROP TABLE IF EXISTS project CASCADE;
CREATE TABLE project (
       db_accession		varchar(48) PRIMARY KEY,
       db_project_link		varchar(80),
       title			text,
       overall_design		text,
       summary			text
);

DROP TABLE IF EXISTS publication CASCADE;
CREATE TABLE publication (
       pubmed_id		text, 
       db_accession		varchar(48) REFERENCES project (db_accession),
       pubmed_link		text,
       doi			text,
       title			text      
);


DROP TABLE IF EXISTS platform CASCADE;
CREATE TABLE platform (
       geo_platform_id		varchar(20) PRIMARY KEY,
       platform_link		varchar(100),
       platform_name		varchar(100),
       annotation_platform  	varchar (30)
);

DROP TABLE IF EXISTS platform_project CASCADE;
CREATE TABLE platform_project (
       geo_platform_id		varchar(20) REFERENCES platform (geo_platform_id),
       db_accession		varchar(48) REFERENCES project (db_accession),
       PRIMARY KEY (geo_platform_id, db_accession)
);

DROP TABLE IF EXISTS transcript CASCADE;
CREATE TABLE transcript (
       transcript_id    varchar(50) PRIMARY KEY     
);


DROP TABLE IF EXISTS platform_transcript CASCADE;
CREATE TABLE platform_transcript (
       geo_platform_id  varchar(20) REFERENCES platform (geo_platform_id),
       transcript_id    varchar(50) REFERENCES transcript (transcript_id)   

);

DROP TABLE IF EXISTS gene_entrez CASCADE;
CREATE TABLE gene_entrez (
       gene_entrez_id    numeric PRIMARY KEY,       
       --gene_symbol	varchar (30),
       species          varchar (50),
       entrez_id_link	text

);

DROP TABLE IF EXISTS transcript_gene_entrez CASCADE;
CREATE TABLE transcript_gene_entrez (
       transcript_id     varchar(50) REFERENCES transcript (transcript_id),            
       gene_entrez_id    numeric REFERENCES gene_entrez (gene_entrez_id)
);


DROP TABLE IF EXISTS gene_symbol CASCADE;
CREATE TABLE gene_symbol (
       gene_symbol_id   varchar(50) PRIMARY KEY,
--       gene_name_desc	text,
       genecards_link	text            
);



DROP TABLE IF EXISTS entrez_symbol CASCADE;
CREATE TABLE entrez_symbol (
       gene_entrez_id    numeric REFERENCES gene_entrez (gene_entrez_id),
       gene_symbol_id    varchar(50) REFERENCES gene_symbol (gene_symbol_id)            
);

DROP TABLE IF EXISTS experiment CASCADE;
CREATE TABLE experiment (
       experiment_id		varchar(80) PRIMARY KEY,
       db_accession		varchar(48) REFERENCES project (db_accession),
       n_samples		numeric,
       sample_description	text,
       observation		text 
);

DROP TABLE IF EXISTS deg_analysis CASCADE;
CREATE TABLE deg_analysis (
       experiment_id               varchar(80) REFERENCES experiment (experiment_id),
       transcript_id               varchar(50) REFERENCES transcript (transcript_id),
       entrez_id                   varchar(80),
       hgnc_symbol                 varchar(200),
       log_fc                      numeric,
       adj_pvalue                  numeric
);

DROP TABLE IF EXISTS experimental_factor_reference CASCADE;
CREATE TABLE experimental_factor_reference  (
       experiment_id               		varchar(48) REFERENCES experiment (experiment_id),
       reference_factor_name              	varchar(80),
       reference_factor_ontology_name     	varchar(80),
       reference_factor_ontology_id       	varchar(80),
       reference_factor_ontology_link     	varchar(200),
       reference_factor_ontology_description	text,
       reference_sp                       	varchar(80),
       reference_sp_id                    	varchar(80),
       reference_tissue                   	varchar(80),
       reference_tissue_ontology          	varchar(80),
       reference_virus                    	varchar(80),
       reference_samples                  	varchar(80),
       PRIMARY KEY (experiment_id, reference_factor_name)
);

DROP TABLE IF EXISTS experimental_factor_test CASCADE;
CREATE TABLE experimental_factor_test  (
       experiment_id               	  varchar(48) REFERENCES experiment (experiment_id),
       test_factor_name                   varchar(80),
       test_factor_ontology_name          varchar(80),
       test_factor_ontology_id            varchar(80),
       test_factor_ontology_link          varchar(200),
       test_factor_ontology_description	  text,
       test_sp                            varchar(80),
       test_sp_id                         varchar(80),
       test_tissue                        varchar(80),
       test_tissue_ontology               varchar(80),
       test_virus                         varchar(80),
       test_samples                       varchar(80),
       PRIMARY KEY (experiment_id, test_factor_name)
);


DROP TABLE IF EXISTS go_analysis CASCADE;
CREATE TABLE go_analysis (
       experiment_id        varchar(80) REFERENCES experiment (experiment_id),
       go_id                varchar(80),
       go_term	       	    text,
       p_value              varchar(80),
       no_go_size           varchar(80),
       entrez		    text,
       PRIMARY KEY (experiment_id, go_id)
);


DROP TABLE IF EXISTS orthologs CASCADE;
CREATE TABLE orthologs(
       entrez_id_mmulatta 		numeric REFERENCES gene_entrez (gene_entrez_id),
       entrez_id_hsapiens 		numeric REFERENCES gene_entrez (gene_entrez_id),
       PRIMARY KEY (entrez_id_mmulatta, entrez_id_hsapiens)
);

DROP TABLE IF EXISTS ontology CASCADE;
CREATE TABLE ontology(
       db_accession 		varchar(48) REFERENCES project (db_accession),
       term			text,
       ontology_name		varchar(50),
       ontology_id		varchar(50),
       tag			varchar(30),   
       PRIMARY KEY (db_accession, term)
);

-- ===============
-- ALTER TABLES --
-- ===============

ALTER TABLE go_analysis ADD FOREIGN KEY (experiment_id) REFERENCES experiment (experiment_id);

ALTER TABLE experimental_factor_reference ADD FOREIGN KEY (experiment_id) REFERENCES experiment (experiment_id);

ALTER TABLE experimental_factor_test ADD FOREIGN KEY (experiment_id) REFERENCES experiment (experiment_id);

ALTER TABLE orthologs ADD FOREIGN KEY (entrez_id_mmulatta) REFERENCES gene_entrez (gene_entrez_id);

ALTER TABLE orthologs ADD FOREIGN KEY (entrez_id_hsapiens) REFERENCES gene_entrez (gene_entrez_id);

-- ========
-- INDEX --
-- ========

--CREATE INDEX ON  gene_symbol (gene_symbol);

CREATE INDEX ON  deg_analysis (log_fc);

CREATE INDEX ON  deg_analysis (adj_pvalue);

CREATE INDEX ON  go_analysis (go_term);

-- ================
-- COPY * TABLES --
-- ================

\copy project FROM '/var/tmp/hihisiv_init/project.csv' DELIMITER ',' CSV

\copy publication FROM '/var/tmp/hihisiv_init/publication.csv' DELIMITER ',' CSV

\copy platform FROM '/var/tmp/hihisiv_init/platform.csv' DELIMITER ',' CSV 

\copy platform_project FROM '/var/tmp/hihisiv_init/platform_project.csv' DELIMITER ',' CSV

\copy transcript FROM '/var/tmp/hihisiv_init/transcript.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy platform_transcript FROM '/var/tmp/hihisiv_init/platform_transcript.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy gene_entrez FROM '/var/tmp/hihisiv_init/gene_entrez.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy gene_symbol FROM '/var/tmp/hihisiv_init/gene_symbol.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy experiment FROM '/var/tmp/hihisiv_init/experiment.csv' DELIMITER ',' CSV

\copy experimental_factor_reference FROM '/var/tmp/hihisiv_init/experimental_factor_reference.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy experimental_factor_test FROM '/var/tmp/hihisiv_init/experimental_factor_test.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy entrez_symbol FROM '/var/tmp/hihisiv_init/entrez_symbol.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy transcript_gene_entrez FROM '/var/tmp/hihisiv_init/transcript_gene_entrez.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy gene_entrez FROM '/var/tmp/hihisiv_init/additional_gene_entrez.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy orthologs FROM '/var/tmp/hihisiv_init/orthologs.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy ontology FROM '/var/tmp/hihisiv_init/ontology.csv' WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy deg_analysis FROM PROGRAM 'awk FNR-1 /var/tmp/hihisiv_init/deg_analysis/*.csv'  WITH (FORMAT CSV, HEADER, DELIMITER ',')

\copy go_analysis FROM PROGRAM 'awk FNR-1 /var/tmp/hihisiv_init/go_analysis/*.csv'  WITH (FORMAT CSV, HEADER, DELIMITER ',')


-- ====================
-- MATERIALIZED VIEW --
-- ====================


