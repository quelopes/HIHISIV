import streamlit as st
import streamlit.components.v1 as components
import psycopg2
import pandas as pd
import networkx as nx
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

#@st.cache
def convert_df(df):
    return df.to_csv().encode('utf-8')

#st.sidebar.image('HIHISIV_logo.png')

st.sidebar.markdown("[![HIHISIV Logo](https://hihisiv.github.io/assets/img/title.png)](https://hihisiv.github.io)", unsafe_allow_html=False)

# Initialize connection.
# Uses st.cache to only run once.
#@st.cache(allow_output_mutation=True, hash_funcs={"_thread.RLock": lambda _: None})
def init_connection():
    return psycopg2.connect(**st.secrets["postgres"])

conn = init_connection()
cur = conn.cursor()

cur.execute("select round(min(logFC),2), round(max(logFC),2) from analysis;")
ranges = cur.fetchall()
minlogfc=float(ranges[0][0])
maxlogfc=float(ranges[0][1])

query_type = st.sidebar.selectbox('Query type:', ('Genes', 'Gene Ontology', 'Transcripts', 'Traits', 'Single gene co-expression network', 'Gene set co-expression network'))
adjPvalue = st.sidebar.radio('Adjusted p-value:', ('0.01','0.05','1.0'), format_func=format_value)
logFC = st.sidebar.slider('Log-FC:', minlogfc, maxlogfc, 0.0)
if(query_type != 'Single gene co-expression network' and query_type != 'Gene set co-expression network'):
    host_type = st.sidebar.radio('Host type:', ('any','non-natural','natural', 'non-natural_vs_natural','human'))

# =============
# === GENES ===
# =============
if query_type=='Genes':
    #cur.execute("select * from gene_symbols;")

    id_type = st.sidebar.radio('Query gene by:', ('Gene symbol', 'Entrez Gene ID'))

    if(id_type == 'Gene symbol'):
        gene_sel = 'gene_symbol'
    else:
        gene_sel = 'entrez_id'

    sql =  "select distinct(" + gene_sel + ") from experiment_gene "
    sql += "where adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s order by " + gene_sel
    cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC})
    gene_sels = cur.fetchall()
    gene_sels_df = pd.DataFrame(gene_sels)
    gene_sel_id = st.sidebar.selectbox(id_type, gene_sels_df)

    # Get corresponding gene_symbol for naming the downloaded csv file
    #cur.execute("select gene_symbol from gene where entrez_id=%(entrez_id)s", {'entrez_id': entrez_id})
    #genes = cur.fetchall()
    #gene_symbol = genes[0][0]

    sql =  "select distinct(tissue) from experiment_gene "
    sql += "where adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s and "
    sql += gene_sel + "=%(gene_sel)s order by tissue"
    cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC, 'gene_sel': gene_sel_id})
    tissues = cur.fetchall()
    tissue_df = pd.DataFrame(tissues)
    tissue_df.loc[-1] = ['any']
    tissue_df.index = tissue_df.index + 1
    tissue_df = tissue_df.sort_index()
    tissue = st.sidebar.selectbox('Tissue:', tissue_df)

    sql =  "select entrez_id, gene_symbol, transcript_id, experiment_id, design_type, tissue, host_type, "
    sql += "to_char(adjPvalue, '9.99EEEE') as adjPvalue, to_char(logFC, '999D99') as logFC "
    sql += "from experiment_gene "
    sql += "where adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s and "
    sql += gene_sel + "=%(gene_sel)s"
    if(host_type != 'any'):
        sql += " and host_type=%(host_type)s"
    if(tissue != 'any'):
        sql += " and tissue=%(tissue)s"

    cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC, 'gene_sel': gene_sel_id, 'tissue': tissue, 'host_type': host_type})
    colnames = [desc[0] for desc in cur.description]
    rows = cur.fetchall()
    df = pd.DataFrame(rows, columns=colnames)
    st.dataframe(df)


    csv = convert_df(df)
    st.download_button(
        label = "Download data as CSV",
        data = csv,
        file_name = str(gene_sel_id) + '.csv',
        mime = 'text/csv',
    )

# =============
# === GENES ===
# =============
if query_type=='Transcripts':
    #cur.execute("select * from gene_symbols;")

    #id_type = st.sidebar.radio('Query transcript by:', ('Gene symbol', 'Entrez Gene ID'))

    sql =  "select distinct(transcript_id) from experiment_gene "
    sql += "where adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s order by transcript_id"
    cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC})
    transcript_sels = cur.fetchall()
    transcript_sels_df = pd.DataFrame(transcript_sels)
    transcript_sel_id = st.sidebar.selectbox("Transcript Id.:", transcript_sels_df)

    # Get corresponding gene_symbol for naming the downloaded csv file
    #cur.execute("select gene_symbol from gene where entrez_id=%(entrez_id)s", {'entrez_id': entrez_id})
    #genes = cur.fetchall()
    #gene_symbol = genes[0][0]

    sql =  "select distinct(tissue) from experiment_gene "
    sql += "where adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s and "
    sql += "transcript_id=%(transcript_sel)s order by tissue"
    cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC, 'transcript_sel': transcript_sel_id})
    tissues = cur.fetchall()
    tissue_df = pd.DataFrame(tissues)
    tissue_df.loc[-1] = ['any']
    tissue_df.index = tissue_df.index + 1
    tissue_df = tissue_df.sort_index()
    tissue = st.sidebar.selectbox('Tissue:', tissue_df)

    sql =  "select entrez_id, gene_symbol, transcript_id, experiment_id, design_type, tissue, host_type, "
    sql += "to_char(adjPvalue, '9.99EEEE') as adjPvalue, to_char(logFC, '999D99') as logFC "
    sql += "from experiment_gene "
    sql += "where adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s and "
    sql += "transcript_id=%(transcript_sel)s"
    if(host_type != 'any'):
        sql += " and host_type=%(host_type)s"
    if(tissue != 'any'):
        sql += " and tissue=%(tissue)s"

    cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC, 'transcript_sel': transcript_sel_id, 'tissue': tissue, 'host_type': host_type})
    colnames = [desc[0] for desc in cur.description]
    rows = cur.fetchall()
    df = pd.DataFrame(rows, columns=colnames)
    st.dataframe(df)


    csv = convert_df(df)
    st.download_button(
        label = "Download data as CSV",
        data = csv,
        file_name = str(transcript_sel_id) + '.csv',
        mime = 'text/csv',
    )

# =====================
# === GENE ONTOLOGY ===
# =====================
if query_type=='Gene Ontology':
    go_domain = st.sidebar.radio('Go domain: ', ('biological_process', 'molecular_function', 'cellular_component'))
    go_type = st.sidebar.radio('Query by:', ('GO term', 'GO ID'))
    if(go_type == 'GO term'):
        sql =  "select distinct(go_term) "
        sql += "from   experiment natural join analysis natural join transcript_gene "
        sql += "natural join gene natural join gene_go natural join go "
        sql += "where adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s and go_domain=%(go_domain)s;"
        cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC, 'go_domain': go_domain})
        go_terms = cur.fetchall()
        go_terms_df = pd.DataFrame(go_terms)
        go_term = st.sidebar.selectbox('GO term:', go_terms_df)
        cur.execute("select go_id from go where go_domain=%(go_domain)s and go_term=%(go_term)s;", {'go_domain': go_domain, 'go_term': go_term})
        go_ids = cur.fetchall()
        go_id = go_ids[0][0]
    else:
        sql =  "select distinct(go_id) "
        sql += "from   experiment natural join analysis natural join transcript_gene "
        sql += "natural join gene natural join gene_go natural join go "
        sql += "where adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s and go_domain=%(go_domain)s;"
        cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC, 'go_domain': go_domain})
        go_ids = cur.fetchall()
        go_ids_df = pd.DataFrame(go_ids)
        go_id = st.sidebar.selectbox('GO ID:', go_ids_df)


    sql =  "select distinct experiment_id, gene_symbol, entrez_id, transcript_id, to_char(adjPvalue, '9.99EEEE') as adjPvalue, to_char(logFC, '999D99') as logFC "
    sql += "from   experiment natural join analysis natural join transcript_gene natural join gene natural join gene_go natural join go where "
    if(host_type != 'any'):
        sql += "host_type=%(host_type)s and "
    if(go_type == 'GO term'):
        sql += "go_term=%(go_term)s and go_domain=%(go_domain)s and adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s;"
        cur.execute(sql, {'go_term': go_term, 'go_domain': go_domain, 'adjPvalue': adjPvalue, 'logFC': logFC, 'host_type': host_type})
    else:
        sql += "go_id=%(go_id)s and go_domain=%(go_domain)s and adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s;"
        cur.execute(sql, {'go_id': go_id, 'go_domain': go_domain, 'adjPvalue': adjPvalue, 'logFC': logFC, 'host_type': host_type})

    colnames = [desc[0] for desc in cur.description]
    rows = cur.fetchall()
    df = pd.DataFrame(rows, columns=colnames)
    st.dataframe(df)

    csv = convert_df(df)
    st.download_button(
        label="Download data as CSV",
        data=csv,
        file_name= go_id + '.csv',
        mime='text/csv',
    )

# ==============
# === TRAITS ===
# ==============
if query_type=='Traits':
    sql =  "select distinct(trait_id) "
    sql += "from   gene_trait natural join gene natural join transcript_gene natural join "
    sql += "analysis "
    sql += "where adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s;"
    cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC})
    trait_ids = cur.fetchall()
    trait_ids_df = pd.DataFrame(trait_ids)
    trait_id = st.sidebar.selectbox('Trait id.:', trait_ids_df)

    sql =  "select distinct(design_type) "
    sql += "from   gene_trait natural join gene natural join transcript_gene natural join "
    sql += "analysis natural join experiment "
    sql += "where adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s;"
    cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC})
    design_types = cur.fetchall()
    design_types_df = pd.DataFrame(design_types)
    design_types_df.loc[-1] = ['any']
    design_types_df.index = design_types_df.index + 1
    design_types_df = design_types_df.sort_index()

    design_type = st.sidebar.selectbox('Design_types:', design_types_df)

    sql =  "select gene_symbol, entrez_id, transcript_id, experiment_id, design_type, host_type, to_char(adjPvalue, '9.99EEEE'), to_char(logFC, '999D99') as logFC "
    sql += "from   gene_trait natural join gene natural join  transcript_gene natural join analysis natural join experiment "
    sql += "where  trait_id=%(trait_id)s and adjPvalue<=%(adjPvalue)s and logFC>=%(logFC)s"
    if(design_type != 'any'):
        sql += " and design_type=%(design_type)s"
    if(host_type != 'any'):
        sql += " and host_type=%(host_type)s"


    cur.execute(sql, {'trait_id': trait_id, 'adjPvalue': adjPvalue, 'logFC': logFC, 'design_type': design_type, 'host_type': host_type})
    colnames = [desc[0] for desc in cur.description]
    rows = cur.fetchall()
    df = pd.DataFrame(rows, columns=colnames)
    st.dataframe(df)

    csv = convert_df(df)
    st.download_button(
        label="Download data as CSV",
        data=csv,
        file_name = trait_id + 'csv',
        mime='text/csv',
    )

# ==============================
# === CO-EXPRESSION ONE GENE ===
# ==============================
if query_type=='Single gene co-expression network':
    sql =  "select distinct gene_symbol from experiment_gene "
    sql += "where adjPvalue<=%(adjPvalue)s and abs(logFC)>=%(logFC)s order by gene_symbol;"
    cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC})
    genes = cur.fetchall()
    genes_df = pd.DataFrame(genes)
    gene_symbol = st.sidebar.selectbox('Gene symbol:', genes_df)
    #ppi = st.sidebar.checkbox('PPI from String db')

    sql =  "select distinct experiment_id, gene_symbol, entrez_id, transcript_id, "
    sql += "to_char(logFC, '999D99') as logFC, to_char(adjPvalue, '9.99EEEE') as adjPvalue "
    sql += "from   experiment_gene "
    sql += "where  adjPvalue<=%(adjPvalue)s and abs(logFC)>=1 and gene_symbol!=%(gene_symbol)s and experiment_id in "
    sql += "(select experiment_id "
    sql += " from   experiment_gene "
    sql += " where  gene_symbol=%(gene_symbol)s and adjPvalue<=%(adjPvalue)s and abs(logFC)>=%(logFC)s) "
    sql += "order by experiment_id;"
    cur.execute(sql, {'gene_symbol': gene_symbol, 'adjPvalue': adjPvalue, 'logFC': logFC})
    colnames = [desc[0] for desc in cur.description]
    rows = cur.fetchall()
    df = pd.DataFrame(rows, columns=colnames)

    sql =  "select gene_symbol, count(*) as value from experiment_gene "
    sql += "where  adjPvalue<=%(adjPvalue)s and abs(logFC)>=1 and "
    sql += "gene_symbol!=%(gene_symbol)s and gene_symbol!='NA' and experiment_id in "
    sql += "(select experiment_id "
    sql += " from   experiment_gene "
    sql += " where  gene_symbol=%(gene_symbol)s and adjPvalue<=%(adjPvalue)s and abs(logFC)>=%(logFC)s) "
    sql += "group by gene_symbol order by value desc limit 50"
    cur.execute(sql, {'gene_symbol': gene_symbol, 'adjPvalue': adjPvalue, 'logFC': logFC})
    colnames = [desc[0] for desc in cur.description]
    rows = cur.fetchall()
    df_net = pd.DataFrame(rows, columns=colnames)

    co_expr_net = Network()

    df_net['selected_gene_symbol'] = gene_symbol

    G = nx.from_pandas_edgelist(df_net, 'selected_gene_symbol', 'gene_symbol', 'value')
    co_expr_net.from_nx(G)

    co_expr_net.show('co_expr_net.html')
    HtmlFile = open("co_expr_net.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    components.html(source_code, height = 700, width=1000)

    st.dataframe(df)

    csv = convert_df(df)
    st.download_button(
        label="Download data as CSV",
        data=csv,
        file_name='co-expression_genes.csv',
        mime='text/csv',
    )

# ===============================
# === CO-EXPRESSION SET GENES ===
# ===============================
if query_type=='Gene set co-expression network':
    sql =  "select distinct entrez_id from experiment_gene "
    sql += "where adjPvalue<=%(adjPvalue)s and abs(logFC)>=%(logFC)s order by entrez_id;"
    cur.execute(sql, {'adjPvalue': adjPvalue, 'logFC': logFC})
    co_threshold = st.sidebar.number_input('Co-ocurrence threshold:', min_value=1, value=1)
    entrez_ids = cur.fetchall()
    entrez_ids_df = pd.DataFrame(entrez_ids)
    entrez_id_sel = st.sidebar.multiselect('Entrez ID:', entrez_ids_df)

    edges = []
    edges_with_colors = []
    colnames = []
    for entrez_id in entrez_id_sel:
        sql =  "select * from (select entrez_id as entrez_id_b, count(*) as value from experiment_gene "
        sql += "where  adjPvalue<=%(adjPvalue)s and abs(logFC)>=1 and "
        sql += "entrez_id!=%(entrez_id)s and experiment_id in "
        sql += "(select experiment_id "
        sql += " from   experiment_gene "
        sql += " where  entrez_id=%(entrez_id)s and adjPvalue<=%(adjPvalue)s and abs(logFC)>=%(logFC)s) "
        sql += "group by entrez_id order by value desc) as s where value>=%(co_threshold)s"
        cur.execute(sql, {'entrez_id': entrez_id, 'adjPvalue': adjPvalue, 'logFC': logFC, 'co_threshold': co_threshold})
        colnames = [desc[0] for desc in cur.description]
        rows = cur.fetchall()
        rows_e = [(entrez_id,) + x for x in filter(lambda tuple: tuple[0] in entrez_id_sel, rows)]
        edges.extend(filter(lambda tuple: tuple[0] < tuple[1], rows_e))
        for tuple in edges:
            sql =  "select entrez_id_a, entrez_id_b from gene_gene "
            sql += "where entrez_id_a=%(node_a)s and entrez_id_b=%(node_b)s or "
            sql += "entrez_id_a=%(node_b)s and entrez_id_b=%(node_a)s"
            cur.execute(sql, {'node_a': tuple[0], 'node_b': tuple[1]})
            if cur.rowcount > 0:
                edges_with_colors.append([tuple[0], tuple[1], tuple[2], 'red'])
            else:
                edges_with_colors.append([tuple[0], tuple[1], tuple[2], 'grey'])

    if edges_with_colors:
        df = pd.DataFrame(edges_with_colors, columns=['entrez_id_a', 'entrez_id_b', 'value', 'color'])
        df_string_id = df.astype({'entrez_id_a': str, 'entrez_id_b': str})

        co_expr_net_mult = Network()
        #G_mult = nx.from_pandas_edgelist(df_string_id, 'entrez_id_a', 'entrez_id_b', 'value')
        G_mult = nx.from_pandas_edgelist(df_string_id, 'entrez_id_a', 'entrez_id_b', ['value', 'color'])
        co_expr_net_mult.from_nx(G_mult)

        co_expr_net_mult.show('co_expr_net_mult.html')
        HtmlFile = open("co_expr_net_mult.html", 'r', encoding='utf-8')
        source_code = HtmlFile.read()
        components.html(source_code, height = 700,width=1000)
        st.text("Red edges indicate that the genes also have a PPI relationship.")
        st.dataframe(df_string_id)


    # sql =  "select distinct entrez_id_a as entrez_id from gene_gene "
    # sql += "union "
    # sql += "select distinct entrez_id_b as entrez_id from gene_gene "
    # sql += "order by entrez_id;"
    # cur.execute(sql)
    # entrez_ids = cur.fetchall()
    # entrez_ids_df = pd.DataFrame(entrez_ids)
    # entrez_id_sel = st.sidebar.multiselect('Entrez ID:', entrez_ids_df)
    #
    # sql =  "select entrez_id_a, entrez_id_b from gene_gene "
    # sql += "where  entrez_id_a=ANY(%(entrez_id_sel)s) and entrez_id_b=ANY(%(entrez_id_sel)s);"
    # cur.execute(sql, {'entrez_id_sel': entrez_id_sel})
    # colnames = [desc[0] for desc in cur.description]
    # rows = cur.fetchall()
    # df = pd.DataFrame(rows, columns=colnames)
    # st.dataframe(df)

        csv = convert_df(df)
        st.download_button(
            label="Download data as CSV",
            data=csv,
            file_name='co-expression_genes.csv',
            mime='text/csv',
        )

    else:
        st.table(edges)
    #co_expr_net = Network()

    #df['selected_gene_symbol'] = gene_symbol

    #G = nx.from_pandas_edgelist(df, 'selected_gene_symbol', 'gene_symbol', 'value')
    #co_expr_net.from_nx(G)

    #co_expr_net.show('co_expr_net.html')
    #HtmlFile = open("co_expr_net.html", 'r', encoding='utf-8')
    #source_code = HtmlFile.read()
    #components.html(source_code, height = 1200,width=1000)
