create table project (
       db_accession		varchar(48) primary key,
       db_project_link	varchar(80),
       bioProject_id		varchar(15),
       title			text,
       overall_design		text,
       summary		text
);

create table publication (
       pubmed_id		text primary key,
       pubmed_link		text,
       doi			text,
       title			text
);

create table platform (
       geo_platform_id		text primary key,
       platform_link		varchar(100),
       platform_name		text
);

create table design (
       design_type		text primary key
);

create table species (
       ncbi_taxonomy_id	numeric primary key,
       taxonomy_id_link	text,
       species_name		varchar(50)
);

create table tissue (
       tissue_name		text primary key
);

create table virus (
       virus_name		text primary key
);

create table experiment (
       experiment_id		varchar(80) primary key,
       db_accession		varchar(48) references project (db_accession),
       design_type		text references design (design_type),
       tissue			text references tissue (tissue_name),
       n_samples		numeric,
       sample_description	text,
       host_type		varchar(25) CHECK (host_type='human' OR host_type='natural' OR host_type='non-natural' OR host_type='non-natural_vs_natural'),
       observation		text
);

create table transcript (
       transcript_id 		varchar(50) primary key
);

create table gene (
       entrez_id		numeric primary key,
       gene_symbol		varchar (30),
       link_ncbi		text,
       gene_name_desc		text,
       chromosome		varchar (50)
);

create table trait (
       trait_id		varchar(30) primary key,
       source_db		text,
       trait_description	text
);

create table go (
       go_id		varchar(30) primary key,
       go_term		text,
       go_domain	varchar(30)
);


create table project_publication (
       pubmed_id		text references publication (pubmed_id),
       db_accession		varchar(48) references project (db_accession),
       primary key (pubmed_id, db_accession)
);

create table experiment_host (
       experiment_id 		varchar(80) references experiment (experiment_id),
       ncbi_taxonomy_id 	numeric references species (ncbi_taxonomy_id),
       gender 			varchar(10) CHECK (gender='male' OR gender='female' OR gender='both' OR gender='undefined' OR gender='NA'),
       age 			text,
       primary key (experiment_id,ncbi_taxonomy_id)
);

create table experiment_virus (
       experiment_id		varchar(80) references experiment (experiment_id),
       virus_name		text references virus (virus_name),
       primary key (experiment_id, virus_name)
);

create table experiment_platform (
       experiment_id		varchar(80) references experiment (experiment_id),
       geo_platform_id		text references platform (geo_platform_id),
       primary key (experiment_id, geo_platform_id)
);

create table platform_transcript (
       platform_id 		text references platform (geo_platform_id),
       transcript_id 		varchar(50) references transcript (transcript_id),
       primary key (platform_id, transcript_id)
);

create table transcript_gene (
       transcript_id 		varchar(50) references transcript (transcript_id),
       entrez_id 		numeric references gene (entrez_id),
       primary key (transcript_id, entrez_id)
);

create table gene_species (
       entrez_id numeric references gene (entrez_id),
       ncbi_taxonomy_id numeric references species (ncbi_taxonomy_id),
       primary key (entrez_id, ncbi_taxonomy_id)
);

create table gene_trait (
        entrez_id		numeric references gene (entrez_id),
	trait_id		varchar(30) references trait (trait_id),
	primary key		(entrez_id, trait_id)
);

create table gene_go (
        entrez_id		numeric references gene (entrez_id),
	go_id  		varchar(30) references go (go_id),
	primary key		(entrez_id, go_id)
);


create table gene_gene (
       entrez_id_a 		numeric references gene (entrez_id),
       entrez_id_b 		numeric references gene (entrez_id),
       relationship 		text,
       source 			text,
       observation text,
       primary key (entrez_id_a, entrez_id_b)
);

create table analysis (
       experiment_id varchar(80) references experiment (experiment_id),
       transcript_id varchar(50) references transcript (transcript_id),
       normalization_method text,
       logFC numeric,
       adjPvalue numeric,
       analysis_script varchar(100),
       primary key (experiment_id, transcript_id)
);

\copy project FROM '/var/tmp/hihisiv_init/project.csv' DELIMITER ',' CSV

\copy publication FROM '/var/tmp/hihisiv_init/publication.csv' DELIMITER ',' CSV

\copy design FROM '/var/tmp/hihisiv_init/design.csv' DELIMITER ',' CSV

\copy species FROM '/var/tmp/hihisiv_init/species.csv' DELIMITER ',' CSV

\copy tissue FROM '/var/tmp/hihisiv_init/tissue.csv' DELIMITER ',' CSV

\copy experiment FROM '/var/tmp/hihisiv_init/experiment.csv' DELIMITER ',' CSV

\copy platform FROM '/var/tmp/hihisiv_init/platform.csv' DELIMITER ',' CSV

\copy virus FROM '/var/tmp/hihisiv_init/virus.csv' DELIMITER ',' CSV

\copy transcript FROM '/var/tmp/hihisiv_init/transcript.csv' DELIMITER ',' CSV

\copy gene FROM '/var/tmp/hihisiv_init/gene.csv' DELIMITER ',' CSV

\copy trait FROM '/var/tmp/hihisiv_init/trait.csv' DELIMITER ',' CSV

\copy go FROM '/var/tmp/hihisiv_init/go.csv' DELIMITER ',' CSV

\copy project_publication FROM '/var/tmp/hihisiv_init/project_publication.csv' DELIMITER ',' CSV

\copy experiment_virus FROM '/var/tmp/hihisiv_init/experiment_virus.csv' DELIMITER ',' CSV

\copy experiment_platform FROM '/var/tmp/hihisiv_init/experiment_platform.csv' DELIMITER ',' CSV

\copy experiment_host FROM '/var/tmp/hihisiv_init/experiment_host.csv' DELIMITER ',' CSV

\copy platform_transcript FROM '/var/tmp/hihisiv_init/platform_transcript.csv' DELIMITER ',' CSV

\copy transcript_gene FROM '/var/tmp/hihisiv_init/transcript_gene.csv' DELIMITER ',' CSV

\copy gene_species FROM '/var/tmp/hihisiv_init/gene_species.csv' DELIMITER ',' CSV

\copy gene_trait FROM '/var/tmp/hihisiv_init/gene_trait.csv' DELIMITER ',' CSV

\copy gene_go FROM '/var/tmp/hihisiv_init/gene_go.csv' DELIMITER ',' CSV

\copy gene_gene FROM '/var/tmp/hihisiv_init/gene_gene.csv' DELIMITER ',' CSV

\copy analysis FROM '/var/tmp/hihisiv_init/analysis.csv' DELIMITER ',' CSV

create materialized view experiment_gene as
  select * from experiment natural join analysis natural join transcript_gene natural join gene where gene_symbol!='';

create materialized view gene_symbols as
  select distinct(gene_symbol) from experiment_gene order by gene_symbol;

create index on analysis (adjPvalue);

create index on analysis (logFC);

create index on gene (gene_symbol);

create index on experiment_gene (logfc);

create index on experiment_gene (adjPValue);

create index on experiment_gene (gene_symbol);

create index on go (go_term);
