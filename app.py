import streamlit as st
import pandas as pd
import os, glob
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple

# Importação dos módulos da arquitetura Enterprise
from modules.processor import VCFProcessor
from modules.visualizer import BioVisualizer
from modules.reporter import ReportManager

# Instanciação dos controladores globais
proc, viz, rep = VCFProcessor(), BioVisualizer(), ReportManager()

def build_agg_table(df: pd.DataFrame, files: List[str]) -> pd.DataFrame:
    """
    Objetivo: Consolidar métricas por amostra para a Tabela Consolidada.
    Lógica: Itera sobre arquivos e agrega contagens e genes únicos por SAMPLEID.
    Saída: pd.DataFrame (risco_maior, genes_encontrados, n_variantes).
    """
    rows = []
    for sid in [f.split('.')[0] for f in files]:
        sub = df[df['SAMPLEID'] == sid]
        rows.append({
            'SAMPLEID': sid,
            'risco_maior': 'SIM' if not sub.empty else 'NÃO',
            'genes_encontrados': ', '.join(sub['GENE'].unique()),
            'n_variantes': len(sub)
        })
    return pd.DataFrame(rows)

def render_tab_geral(df: pd.DataFrame, files: List[str]):
    """
    Objetivo: Renderizar o dashboard principal conforme layout solicitado.
    Lógica: Gráfico de Barras (esquerda), Pizza (direita) e Tabela (abaixo).
    """
    col_l, col_r = st.columns([2, 1])
    with col_l: st.pyplot(viz.plot_gene_bar(df))
    with col_r: st.pyplot(viz.plot_risk_pie(df, files))
    st.divider()
    st.subheader("Tabela Consolidada de Amostras")
    st.dataframe(build_agg_table(df, files), use_container_width=True)

def main():
    """Ponto de entrada: Gerencia Sidebar e Session State para evitar 'tela preta'."""
    st.sidebar.title("Painel de Controle")
    p = { 'dp_min': st.sidebar.slider("DP Min", 10, 200, 20), 'vaf_min': st.sidebar.slider("VAF Min", 0.0, 1.0, 0.05),
          'max_pop_af': st.sidebar.slider("gnomAD Max AF", 0.0, 0.05, 0.01, format="%.3f"),
          'genes': [i.strip() for i in os.getenv('GENES_ALTO_RISCO', '').split(',')],
          'only_p': st.sidebar.checkbox("Apenas ClinVar Pathogenic"), 'impacts': ['HIGH', 'MODERATE'],
          'cons': [i.strip() for i in os.getenv('CONSEQUENCIAS_INTERESSE', '').split(',')] }

    if st.sidebar.button("Iniciar Análise de VCFs"):
        paths = glob.glob(os.path.join(os.getenv('INPUT_DIR'), "**/*.vcf*"), recursive=True)
        with st.spinner("Processando coorte em paralelo..."):
            st.session_state['df'] = pd.DataFrame(proc.run_parallel(paths, p))
            st.session_state['files'] = [os.path.basename(f) for f in paths]

    if 'df' in st.session_state:
        t1, t2, t3, t4, t5 = st.tabs(["Geral", "OncoPrint", "Lollipop", "Assinatura", "PDF"])
        with t1: render_tab_geral(st.session_state['df'], st.session_state['files'])
        with t2: st.pyplot(viz.plot_oncoprint(st.session_state['df']))
        with t3: 
            g = st.selectbox("Gene", st.session_state['df']['GENE'].unique())
            st.pyplot(viz.plot_lollipop(st.session_state['df'], g))
        with t5:
            sid = st.selectbox("Amostra", st.session_state['df']['SAMPLEID'].unique())
            if sid:
                df_s = st.session_state['df'][st.session_state['df']['SAMPLEID'] == sid]
                pdf = rep.create_pdf(sid, df_s, p, [viz.plot_lollipop(df_s, df_s['GENE'].iloc[0])])
                st.download_button(f"Baixar PDF - {sid}", pdf, f"{sid}.pdf")

if __name__ == "__main__":
    main()