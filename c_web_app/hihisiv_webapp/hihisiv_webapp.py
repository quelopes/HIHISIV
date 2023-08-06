import streamlit as st
import streamlit.components.v1 as components
import psycopg2
import pandas as pd
import networkx as nx
from hihisiv_queries import *
from pyvis.network import Network

st.set_page_config(layout="wide")

hide_st_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            header {visibility: hidden;}
            </style>
            """
st.markdown(hide_st_style, unsafe_allow_html=True)

def format_value(value):
    if(value=='1.0'):
        return "any"
    else:
        return value

def convert_df(df):
    return df.to_csv().encode('utf-8')

def save_csv(df, filename):
    csv = convert_df(df)
    st.download_button(
        label = "Download data as CSV",
        data = csv,
        file_name = filename + '.csv',
        mime = 'text/csv',
    )

@st.cache_data
def execute_query(_cursor, query, parameters):
    _cursor.execute(query, parameters)
    colnames = [desc[0] for desc in _cursor.description]
    rows = _cursor.fetchall()
    df = pd.DataFrame(rows, columns=colnames)
    return df

def display_network(df, source, target, edge_attribute, output_filename):
    co_expr_net = Network()
    G = nx.from_pandas_edgelist(df, source, target, edge_attribute)
    co_expr_net.from_nx(G)
    co_expr_net.show(output_filename, notebook=False)
    HtmlFile = open(output_filename, 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    components.html(source_code, height = 700, width=1000)

def display_network_with_isolated(df, source, target, edge_attribute, single_nodes, output_filename):
    co_expr_net = Network()
    G = nx.Graph()
    for node in single_nodes:
        G.add_node(node)
    edges_to_add = [(row['gene_symbol_a'], row['gene_symbol_b'], {'weight': row['value']}) for _, row in df.iterrows()]
    G.add_edges_from(edges_to_add)
    co_expr_net.from_nx(G)
    co_expr_net.show(output_filename, notebook=False)
    HtmlFile = open(output_filename, 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    components.html(source_code, height = 700, width=1000)

st.sidebar.markdown("[![HIHISIV Logo](https://hihisiv.github.io/assets/img/HIHISIV_logo.png)](https://hihisiv.github.io)", unsafe_allow_html=False)

# Initialize connection.
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])

conn = init_connection()
cur = conn.cursor()

# Get maximum absolute value for log_fc
cur.execute(sql_max_abs_log_fc)
ranges = cur.fetchall()
max_abs_log_fc=float(ranges[0][0])

# Selection of query types
query_type = st.sidebar.selectbox('Query type:', ('Genes', 'Transcripts', 'Biological Process (GO)', 'Ontology Terms', 'Single gene co-expression network', 'Gene set co-expression network'))

# Displays adj_pvalue and log_fc parameters only on some query types
if query_type in ['Genes', 'Transcripts', 'Single gene co-expression network', 'Gene set co-expression network']:
#    adjPvalue = st.sidebar.radio('Adjusted p-value:', ('0.01','0.05','1.0'), index=1, format_func=format_value)
    adjPvalue = st.sidebar.slider('Adjusted p-value ≤:', 0.00, 1.00, 0.05, 0.01)
    logFC = st.sidebar.slider('|Log-FC| ≥:', 0.0, max_abs_log_fc, 1.0, 0.1, '%.1f')

# =============
# === GENES ===
# =============
if query_type=='Genes':
    id_type = st.sidebar.radio('Query gene by:', ('Gene symbol', 'Entrez Gene ID'))

    if(id_type == 'Gene symbol'):
        gene_sels_df = execute_query(cur, sql_list_genes_by_symbol, {'adjPvalue': adjPvalue, 'logFC': logFC})
    else:
        gene_sels_df = execute_query(cur, sql_list_genes_by_entrez_id, {'adjPvalue': adjPvalue, 'logFC': logFC})

    gene_sel_id = st.sidebar.selectbox(id_type, gene_sels_df)


    if(id_type == 'Gene symbol'):
        df = execute_query(cur, sql_gene_expression_experiment_by_symbol, 
                    {'adjPvalue': adjPvalue, 'logFC': logFC, 'gene_sel_id': gene_sel_id})
    else:
        df = execute_query(cur, sql_gene_expression_experiment_by_entrez_id, 
                    {'adjPvalue': adjPvalue, 'logFC': logFC, 'gene_sel_id': gene_sel_id})
    df.n_samples = df.n_samples.astype(int)
    st.dataframe(df)

    save_csv(df, str(gene_sel_id))

# ===================
# === TRANSCRIPTS ===
# ===================

if query_type=='Transcripts':
    transcript_sels_df = execute_query(cur, sql_list_transcripts, {'adjPvalue': adjPvalue, 'logFC': logFC})
    transcript_sel_id = st.sidebar.selectbox("Transcript Id.", transcript_sels_df)

    df = execute_query(cur, sql_transcript_experiment, 
                    {'adjPvalue': adjPvalue, 'logFC': logFC, 'transcript_sel_id': transcript_sel_id})

    df.n_samples = df.n_samples.astype(int)
    st.dataframe(df)

    save_csv(df, str(transcript_sel_id))

# ==========================
# === BIOLOGICAL PROCESS ===
# ==========================
if query_type=='Biological Process (GO)':
    go_terms_df = execute_query(cur, sql_list_go_terms, {})
    go_term = st.sidebar.selectbox('GO term:', go_terms_df)

    cur.execute(sql_get_go_id, {'go_term': go_term})
    go_ids = cur.fetchall()
    go_id = go_ids[0][0]

    df = execute_query(cur, sql_experiment_go_term, {'go_term': go_term})

    st.dataframe(df)

    save_csv(df, go_id)

# ==============================
# === CO-EXPRESSION ONE GENE ===
# ==============================

if query_type=='Single gene co-expression network':
    # select distinct relevant gene symbols or ids for use in the gene selection menu
    genes_df = execute_query(cur, sql_list_genes_by_symbol, {'adjPvalue': adjPvalue, 'logFC': logFC})
    gene_symbol = st.sidebar.selectbox('Gene symbol:', genes_df)

    df_net = execute_query(cur, sql_single_gene_co_expression_graph, 
                    {'gene_symbol': gene_symbol, 'adjPvalue': adjPvalue, 'logFC': logFC})

    df_net['selected_gene_symbol'] = gene_symbol

    display_network(df_net, 'selected_gene_symbol', 'hgnc_symbol', 'value', 'co_expr_net.html')
    st.text("Edge thickness in the graph indicates the frequency in which the genes appear co-expressed.")
    df_net = df_net.rename(columns={"value": "frequency"})

    st.dataframe(df_net[['selected_gene_symbol', 'hgnc_symbol', 'experiment_list', 'frequency']])

    save_csv(df_net, 'co-expression_genes')

#===============================
#=== CO-EXPRESSION SET GENES ===
#===============================
if query_type=='Gene set co-expression network':
    display_isolated = st.sidebar.radio('Display isolated nodes:', ('Yes', 'No'))
    genes_df = execute_query(cur, sql_list_genes_by_symbol, {'adjPvalue': adjPvalue, 'logFC': logFC})
    co_threshold = st.sidebar.number_input('Co-ocurrence threshold:', min_value=1, value=1)
    genes_sel = st.sidebar.multiselect('Gene symbol:', genes_df)

    edges = []
    edges_with_colors = []
    colnames = []
    for gene_symbol in genes_sel:
        cur.execute(sql_co_expressed_genes, 
                    {'gene_symbol': gene_symbol, 'adjPvalue': adjPvalue, 'logFC': logFC, 'co_threshold': co_threshold})
        colnames = [desc[0] for desc in cur.description]
        rows = cur.fetchall()
        rows_e = [(gene_symbol,) + x for x in filter(lambda tuple: tuple[0] in genes_sel, rows)]
        edges.extend(filter(lambda tuple: tuple[0] < tuple[1], rows_e))

    df = pd.DataFrame(edges, columns=['gene_symbol_a', 'gene_symbol_b', 'value', 'experiment_list'])
    df_string_id = df.astype({'gene_symbol_a': str, 'gene_symbol_b': str})
    if display_isolated == 'Yes':
        display_network_with_isolated(df_string_id, 'gene_symbol_a', 'gene_symbol_b', ['value'], genes_sel, 'co_expr_net_set.html')
    # check if edges is non-empty
    elif edges:
        display_network(df_string_id, 'gene_symbol_a', 'gene_symbol_b', ['value'], 'co_expr_net_set.html')
    
    st.text("Edge thickness in the graph indicates the frequency in which the genes appear co-expressed.")
    df = df.rename(columns={"value": "frequency"}).sort_values(by=['frequency', 'gene_symbol_a', 'gene_symbol_b'], ascending=False)
    st.dataframe(df[['gene_symbol_a', 'gene_symbol_b', 'experiment_list', 'frequency']])

    save_csv(df, 'co-expression_gene_set')

# ========================================================
# === MODE 6 experiments related to the ontology terms ===
# ========================================================
if query_type=='Ontology Terms':
    terms_sels_df = execute_query(cur, sql_list_ontology_terms, {})
    term_sel = st.sidebar.selectbox("Ontology term", terms_sels_df)

    df = execute_query(cur, sql_experiment_ontology_terms, {'term_sel': term_sel})

    st.dataframe(df)

    save_csv(df, str(term_sel))
