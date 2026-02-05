import io, pandas as pd  # Importação de bibliotecas para buffer de memória e estruturas de dados tabulares
from fpdf import FPDF  # Importação da biblioteca FPDF para construção de documentos PDF

class ReportManager:
    '''
    Descrição: Classe responsável pela exportação do Dossiê Mestre da Coorte com interpretações detalhadas.
    Lógica: Centraliza gráficos e insere blocos de texto que explicam a ciência por trás de cada visualização e seus resultados.
    '''

    def _convert_fig_to_buffer(self, fig):
        '''
        Descrição: Converte um objeto Figure do Matplotlib num buffer de imagem binária PNG.
        Parâmetros:
            - fig (matplotlib.figure.Figure): Instância do gráfico gerado.
        Tipo: Objeto de imagem Matplotlib.
        Entrada: Figure.
        Saída: io.BytesIO contendo a imagem.
        Lógica: Salva a figura no buffer com DPI 150 e reseta o ponteiro de leitura para o início.
        '''
        buffer = io.BytesIO()  # Inicializa o buffer de memória virtual
        fig.savefig(buffer, format='png', bbox_inches='tight', dpi=150)  # Salva a imagem como PNG otimizado
        buffer.seek(0)  # Move o cursor de leitura para o início do buffer
        return buffer  # Retorna o buffer pronto para inserção no PDF

    def _draw_table(self, pdf, df):
        '''
        Descrição: Renderiza a tabela de risco centralizada com explicações sobre os campos.
        Parâmetros:
            - pdf (fpdf.FPDF): Instância do documento PDF.
            - df (pd.DataFrame): Tabela de classificação de risco.
        Tipo: Objeto FPDF e DataFrame Pandas.
        Entrada: Instância do PDF e dados.
        Saída: Nenhuma (modifica o PDF).
        Lógica: Desenha células centralizadas com largura fixa de 40mm para cada coluna técnica.
        '''
        pdf.set_font("Helvetica", 'B', 8)  # Define fonte negrito para o cabeçalho
        col_w, table_w = 40, 160  # Define larguras (4 colunas de 40mm = 160mm total)
        start_x = (210 - table_w) / 2  # Calcula o ponto X para centralização em página A4 (210mm)
        cols = ['SAMPLEID', 'MAIOR_RISCO', 'TP53_PRESENTE', 'N_VARIANTES']  # Colunas selecionadas
        pdf.set_x(start_x)  # Posiciona no X calculado
        for c in cols: pdf.cell(col_w, 8, c, border=1, align='C')  # Desenha células do cabeçalho
        pdf.ln(); pdf.set_font("Helvetica", '', 7)  # Salta linha e altera fonte para os dados
        for _, row in df.iterrows():  # Itera sobre cada amostra processada
            pdf.set_x(start_x)  # Reseta o cursor no X centralizado para cada nova linha
            pdf.cell(col_w, 6, str(row['SAMPLEID']), 1, 0, 'C')  # Célula ID
            pdf.cell(col_w, 6, str(row['MAIOR_RISCO']), 1, 0, 'C')  # Célula Status Risco
            pdf.cell(col_w, 6, str(row['TP53_PRESENTE']), 1, 0, 'C')  # Célula Status TP53
            pdf.cell(col_w, 6, str(row.get('N_VARIANTES', 0)), 1, 1, 'C')  # Célula Contagem final

    def create_cohort_report(self, p: dict, figs: dict, summary: str, df_risk: pd.DataFrame) -> bytes:
        '''
        Descrição: Orquestra a criação do reporter completo com textos explicativos de resultado e centralização.
        Parâmetros:
            - p (dict): Thresholds de filtragem técnica.
            - figs (dict): Dicionário de objetos Figure.
            - summary (str): Texto de resumo executivo.
            - df_risk (pd.DataFrame): Tabela de classificação de risco.
        Tipo: Dict, Dict, String e DataFrame.
        Entrada: Metadados, gráficos e dados.
        Saída: bytes contendo o PDF binário.
        Lógica: Adiciona páginas e insere explicações sobre o significado clínico de cada gráfico antes da imagem centralizada.
        '''
        pdf = FPDF()  # Instancia o motor de PDF
        pdf.set_auto_page_break(auto=True, margin=15)  # Ativa quebra automática de página
        pdf.add_page()  # Adiciona a primeira página (A4)
        
        # --- TÍTULO E RESUMO EXECUTIVO ---
        pdf.set_font("Helvetica", 'B', 16)  # Fonte para o título principal
        pdf.cell(0, 15, "ANÁLISE GENÔMICA MF", ln=True, align='C')  # Título centralizado
        
        pdf.ln(5)  # Espaçamento vertical
        pdf.set_font("Helvetica", 'B', 12); pdf.set_text_color(200, 0, 0)  # Fonte vermelha para destaque
        pdf.cell(0, 10, summary, ln=True, align='C')  # Exibe resumo: Amostras | Risco | TP53
        pdf.set_text_color(0, 0, 0)  # Restaura cor preta
        
        # --- CONFIGURAÇÕES CENTRALIZADAS ---
        pdf.set_font("Helvetica", 'B', 10); pdf.cell(0, 8, "CONFIGURAÇÕES DE FILTRAGEM", ln=True, align='C')  # Título seção
        pdf.set_font("Helvetica", '', 9)  # Fonte metadados
        txt_p = f"DP Mínimo: {p['dp_min']} | VAF Mínimo: {p['vaf_min']} | gnomAD Máximo: {p['max_pop_af']}"  # Texto filtros
        pdf.cell(0, 6, txt_p, ln=True, align='C')  # Thresholds centralizados
        
        # --- EXPLICAÇÃO PIZZA (STATUS DE RISCO) ---
        pdf.ln(5)  # Espaçamento
        pdf.set_font("Helvetica", 'B', 11); pdf.cell(0, 10, "1. STATUS DE RISCO", ln=True)  # Cabeçalho 1
        pdf.set_font("Helvetica", '', 9)  # Fonte descrição
        desc_risk = "O gráfico de pizza representa a penetrância das variantes de alto risco na coorte. Amostras classificadas como 'SIM' possuem pelo menos uma mutação em genes drivers (como JAK2, CALR ou MPL) que atendem aos critérios de patogenicidade."  # Explicação
        pdf.multi_cell(0, 5, desc_risk)  # Texto explicativo
        pdf.image(self._convert_fig_to_buffer(figs['pie']), x=55, w=100)  # Pizza centralizada
        
        # --- EXPLICAÇÃO BARRAS (PREVALÊNCIA) ---
        pdf.add_page()  # Nova página para prevalência
        pdf.set_font("Helvetica", 'B', 11); pdf.cell(0, 10, "2. PREVALÊNCIA POR GENE", ln=True)  # Cabeçalho 2
        pdf.set_font("Helvetica", '', 9)  # Fonte descrição
        desc_prev = "Este gráfico de barras quantifica a frequência absoluta de mutações para cada gene driver. Genes com maior frequência indicam os principais mecanismos de progressão da doença na coorte. O resultado revela quais mutações são dominantes no perfil genômico estudado."  # Explicação
        pdf.multi_cell(0, 5, desc_prev)  # Texto explicativo
        pdf.image(self._convert_fig_to_buffer(figs['gene']), x=25, w=160)  # Barras centralizadas
        
        # --- EXPLICAÇÃO TABELA ---
        pdf.ln(5)  # Espaçamento
        pdf.set_font("Helvetica", 'B', 11); pdf.cell(0, 10, "3. DETALHAMENTO POR AMOSTRA", ln=True)  # Cabeçalho 3
        pdf.set_font("Helvetica", '', 9)  # Fonte descrição
        desc_tab = "A tabela consolidada detalha o status individual. A coluna 'TP53_PRESENTE' é um marcador de extrema relevância clínica; mutações em TP53 são frequentemente associadas à transformação leucêmica e resistência terapêutica em Mielofibrose."  # Explicação
        pdf.multi_cell(0, 5, desc_tab)  # Texto explicativo
        pdf.ln(2); self._draw_table(pdf, df_risk)  # Desenha tabela centralizada
        
        # --- EXPLICAÇÃO ONCOPRINT E ASSINATURAS ---
        pdf.add_page()  # Nova página para análises complexas
        pdf.set_font("Helvetica", 'B', 11); pdf.cell(0, 10, "4. PAISAGEM MUTACIONAL (ONCOPRINT)", ln=True)  # Cabeçalho 4
        pdf.set_font("Helvetica", '', 9)  # Fonte descrição
        desc_onco = "O OncoPrint visualiza a co-ocorrência de mutações entre amostras. O resultado permite identificar se certas mutações ocorrem simultaneamente (co-ocorrência) ou se são mutuamente exclusivas, auxiliando na compreensão da arquitetura clonal da coorte."  # Explicação
        pdf.multi_cell(0, 5, desc_onco)  # Texto explicativo
        pdf.image(self._convert_fig_to_buffer(figs['onco']), x=15, w=180)  # OncoPrint centralizado
        
        pdf.ln(10)  # Espaçamento
        pdf.set_font("Helvetica", 'B', 11); pdf.cell(0, 10, "5. ASSINATURAS DE SUBSTITUIÇÃO SNV", ln=True)  # Cabeçalho 5
        pdf.set_font("Helvetica", '', 9)  # Fonte descrição
        desc_sign = "As assinaturas mutacionais analisam o tipo químico das trocas de bases (ex: C>T). O resultado ajuda a discernir a etiologia das mutações, como processos de envelhecimento celular ou erros em vias específicas de reparo de DNA."  # Explicação
        pdf.multi_cell(0, 5, desc_sign)  # Texto explicativo
        pdf.image(self._convert_fig_to_buffer(figs['sign']), x=35, w=140)  # Assinaturas centralizadas
        
        return bytes(pdf.output())  # Retorna o PDF em formato binário