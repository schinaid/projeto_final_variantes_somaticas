import streamlit as st
import pandas as pd
import os, glob
from typing import List, Dict, Tuple
from modules.processor import VCFProcessor
from modules.visualizer import BioVisualizer
from modules.reporter import ReportManager

# Instanciação dos módulos Enterprise
proc, viz, rep = VCFProcessor(), BioVisualizer(), ReportManager()

def build_agg_table(df: pd.DataFrame, files: List[str]) -> pd.DataFrame:
    """Objetivo: Tabela Consolidada. Lógica: Agregação por Sample ID."""
    rows = []
    for sid in [f.split('.')[0] for f in files]:
        sub = df[df['SAMPLEID'] == sid]
        rows.append({
            'SAMPLEID': sid, 'risco_maior': 'SIM' if not sub.empty else 'NÃO',
            'genes_encontrados': ', '.join(sub['GENE'].unique()), 'n_variantes': len(sub)
        })
    return pd.DataFrame(rows)

def main():
    """Função Principal: Sidebar, Session State e Abas com Explicações."""
    st.sidebar.title("⚙️ Painel de Controle")
    # Captura de parâmetros dinâmicos via Sidebar
    p = { 'dp_min': st.sidebar.slider("DP Mínimo", 10, 200, 20),
          'vaf_min': st.sidebar.slider("VAF Mínimo", 0.0, 1.0, 0.05),
          'max_pop_af': st.sidebar.slider("gnomAD Max AF", 0.0, 0.05, 0.01, format="%.3f"),
          'genes': [i.strip() for i in os.getenv('GENES_ALTO_RISCO', '').split(',')],
          'only_p': st.sidebar.checkbox("Apenas ClinVar Pathogenic") }

    if st.sidebar.button("🚀 Iniciar Análise"):
        paths = glob.glob(os.path.join(os.getenv('INPUT_DIR'), "**/*.vcf*"), recursive=True)
        with st.spinner("Processando coorte..."):
            st.session_state['df'] = pd.DataFrame(proc.run_parallel(paths, p))
            st.session_state['files'] = [os.path.basename(f) for f in paths]

    if 'df' in st.session_state:
        df, f_list = st.session_state['df'], st.session_state['files']
        t1, t2, t3, t4, t5 = st.tabs(["Geral", "OncoPrint", "Lollipop", "Assinaturas", "PDF"])
        
        with t1:
            st.info("**Visão Geral**: Consolida a frequência de mutações e a classificação de risco de toda a coorte.")
            cl, cr = st.columns([2, 1]); cl.pyplot(viz.plot_gene_bar(df)); cr.pyplot(viz.plot_risk_pie(df, f_list))
            st.divider(); st.subheader("📋 Tabela Consolidada de Amostras"); st.dataframe(build_agg_table(df, f_list), use_container_width=True)
        
        with t2:
            st.markdown("### OncoPrint: Paisagem Mutacional")
            st.write("Esta visualização identifica padrões de **co-ocorrência** ou **exclusividade mútua** entre genes em diferentes pacientes.")
            st.pyplot(viz.plot_oncoprint(df))
        
        with t3:
            st.markdown("### Lollipop Plot: Hotspots Proteicos")
            st.write("Identifica onde as mutações se acumulam na estrutura da proteína. A altura da haste indica a Variant Allele Frequency ($$VAF$$).")
            g = st.selectbox("Selecione o Gene", df['GENE'].unique() if not df.empty else [])
            if g: st.pyplot(viz.plot_lollipop(df, g))
            
        with t4:
            st.markdown("### Assinaturas: Perfil de Substituição")
            st.write("Analisa os processos químicos que causaram as mutações (ex: envelhecimento, reparo de DNA).")
            st.pyplot(viz.plot_signatures(df))
            
        with t5:
            st.markdown("### Exportação de PDF Técnico")
            st.info("Gera um PDF formatado com evidências gráficas e metadados de filtragem.")
            sid = st.selectbox("Amostra", df['SAMPLEID'].unique() if not df.empty else [])
            if sid:
                df_s = df[df['SAMPLEID'] == sid]
                pdf_b = rep.create_pdf(sid, df_s, p, viz.plot_lollipop(df_s, df_s['GENE'].iloc[0]))
                st.download_button(f"Baixar PDF - {sid}", pdf_b, f"{sid}.pdf")

if __name__ == "__main__":
    main()