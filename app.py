import streamlit as st # Importação da biblioteca de interface
import pandas as pd # Estruturas de dados
import os, glob # Operações de arquivo
from modules.processor import VCFProcessor # Motor de bioinformática
from modules.visualizer import BioVisualizer # Camada visual
from modules.reporter import ReportManager # Gestor de laudos

class GenomicApp:
    '''
    Descrição: Controlador principal Enterprise para análise de variantes MF.
    Lógica: Orquestra o fluxo de dados, renderiza a UI e aplica explicações científicas.
    '''

    def __init__(self):
        '''Descrição: Inicializa módulos e UI. Lógica: Injeta dependências e define layout.'''
        self.proc, self.viz, self.rep = VCFProcessor(), BioVisualizer(), ReportManager() # Injeção de dependências
        st.set_page_config(page_title="MF Analyzer", layout="wide") # Configura Streamlit

    def _generate_sample_risk_df(self, df: pd.DataFrame, paths: list) -> pd.DataFrame:
        '''Descrição: Gera tabela consolidada. Lógica: Identifica MAIOR_RISCO e TP53 por amostra.'''
        rows = [] # Inicializa acumulador
        for p in paths: # Itera sobre arquivos processados
            sid = os.path.basename(p).split('.')[0] # Extrai ID
            sub = df[df['SAMPLEID'] == sid] # Filtra variantes da amostra específica
            rows.append({ # Consolida dicionário de flags
                'SAMPLEID': sid, 'MAIOR_RISCO': 'SIM' if not sub.empty else 'NÃO',
                'TP53_PRESENTE': 'SIM' if 'TP53' in sub['GENE'].values else 'NÃO',
                'GENES': ', '.join(sub['GENE'].unique()), 'N_VARIANTES': len(sub)
            })
        return pd.DataFrame(rows) # Retorna estrutura tabular consolidada

    def run(self):
        '''Descrição: Loop de execução Streamlit. Lógica: Sidebar, Tabs e Explicações LaTeX.'''
        st.sidebar.title("Painel de controle MF") # Título Sidebar
        p = { 'dp_min': st.sidebar.slider("DP Mínimo", 10, 200, 20), # Slider DP
              'vaf_min': st.sidebar.slider("VAF Mínimo", 0.0, 1.0, 0.05), # Slider VAF
              'max_pop_af': st.sidebar.slider("gnomAD Máximo", 0.0, 0.05, 0.01, format="%.3f") } # Slider gnomAD

        if st.sidebar.button("🚀 Iniciar Processamento"): # Trigger análise
            paths = glob.glob(os.path.join(os.getenv('INPUT_DIR', './inputs'), "**/*.vcf*"), recursive=True) # Busca VCFs
            with st.spinner("Analisando coorte..."): # Feedback visual
                st.session_state['df'] = pd.DataFrame(self.proc.run_parallel(paths, p)) # Executa motor
                st.session_state['files'] = paths # Armazena caminhos
                st.session_state['df'].to_csv('variants_high_risk.tsv', sep='\t', index=False) # Salva TSV 1

        if 'df' in st.session_state: # Verifica se há dados para exibir
            df, files = st.session_state['df'], st.session_state['files'] # Recupera estado
            df_risk = self._generate_sample_risk_df(df, files) # Gera tabela de risco
            df_risk.to_csv('sample_risk.tsv', sep='\t', index=False) # Salva TSV 2
            summary = f"ANÁLISE: {len(df_risk)} Amostras | {len(df_risk[df_risk['MAIOR_RISCO']=='SIM'])} Risco Alto | {len(df_risk[df_risk['TP53_PRESENTE']=='SIM'])} TP53 Mutado" # Resumo
            
            tabs = st.tabs(["Geral", "OncoPrint", "Assinaturas", "Exportar PDF"]) # Cria abas
            figs = {'gene': self.viz.plot_gene_frequency(df), 'onco': self.viz.plot_oncoprint(df), 
                    'sign': self.viz.plot_signatures(df), 'pie': self.viz.plot_risk_pie(df, files)} # Gera gráficos

            with tabs[0]: # ABA GERAL
                st.markdown("### Prevalência e Classificação")
                st.write("A incidência genômica $Freq$ por gene $G_i$ na coorte de tamanho $n$ é calculada por:")
                st.latex(r"Freq(G_i) = \sum_{j=1}^{n} \mathbb{1}_{G_i \in Sample_j}") # Fórmula LaTeX
                st.code(summary) # Exibe resumo textual
                cl, cr = st.columns([2, 1]); cl.pyplot(figs['gene']); cr.pyplot(figs['pie']) # Renderiza gráficos
                st.dataframe(df_risk, use_container_width=True) # Renderiza tabela

            with tabs[1]: # ABA ONCOPRINT
                st.markdown("### OncoPrint: Paisagem Mutacional")
                st.write("Representação da matriz binária de status mutacional $M_{ij}$ para o gene $i$ na amostra $j$:")
                st.latex(r"M_{ij} = \{1 \text{ se variante presente}, 0 \text{ se selvagem}\}") # Fórmula LaTeX
                st.pyplot(figs['onco']) # Renderiza OncoPrint

            with tabs[2]: # ABA ASSINATURAS
                st.markdown("### Assinaturas de Substituição")
                st.write("Cálculo do perfil SNV baseado na Variant Allele Frequency ($$VAF$$):")
                st.latex(r"VAF = \frac{AD_{Alt}}{AD_{Ref} + AD_{Alt}}") # Fórmula LaTeX
                st.pyplot(figs['sign']) # Renderiza Assinaturas

            with tabs[3]: # ABA EXPORTAÇÃO
                st.markdown("### Exportação") # Markdown explicativo
                st.info("O PDF unifica todos os gráficos e a tabela centralizada.") # Caixa de informação
                # FIX: Inclusão do df_risk como argumento para evitar o TypeError
                pdf_data = self.rep.create_cohort_report(p, figs, summary, df_risk) # Chama criação do PDF
                st.download_button("Baixar Relatório", pdf_data, "MF_Report.pdf") # Botão de download
if __name__ == "__main__":
    GenomicApp().run() # Executa aplicação