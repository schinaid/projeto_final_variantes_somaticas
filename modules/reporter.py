import io
from fpdf import FPDF
import pandas as pd

class ReportManager:
    """Objetivo: Gerenciar a exportação de laudos técnicos em PDF."""

    def _fig_to_buf(self, fig) -> io.BytesIO:
        """Objetivo: Converter plot para imagem binária. Saída: Buffer."""
        buf = io.BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
        buf.seek(0)
        return buf

    def _add_filters_info(self, pdf, p: dict):
        """Objetivo: Inserir seção de parâmetros de filtragem. Lógica: Metadados técnicos."""
        pdf.set_fill_color(240, 240, 240)
        pdf.set_font("Helvetica", 'B', 10)
        pdf.cell(0, 10, " PARÂMETROS DE FILTRAGEM (REFERÊNCIA TÉCNICA)", ln=True, fill=True)
        pdf.set_font("Helvetica", '', 9)
        # Detalha DP, VAF e gnomAD aplicados no pipeline
        txt = f"DP Min: {p['dp_min']} | VAF Min: {p['vaf_min']:.2f} | gnomAD Max: {p['max_pop_af']:.3f}"
        pdf.cell(0, 8, txt, ln=True)
        pdf.ln(5)

    def create_pdf(self, sid: str, df: pd.DataFrame, p: dict, fig=None) -> bytes:
        """
        Objetivo: Gerar o PDF completo com gráficos e dados.
        Entrada: ID Amostra, Variantes, Parâmetros, Figura. Saída: bytes.
        """
        pdf = FPDF()
        pdf.add_page()
        pdf.set_font("Helvetica", 'B', 18)
        pdf.cell(0, 15, f"PDF GENÔMICO: {sid}", ln=True, align='C')
        self._add_filters_info(pdf, p)
        # Insere a evidência gráfica (Lollipop) se fornecida
        if fig: pdf.image(self._fig_to_buf(fig), w=170); pdf.ln(80)
        # Lista detalhadamente cada variante encontrada
        for _, r in df.iterrows():
            pdf.set_font("Helvetica", 'B', 11); pdf.cell(0, 8, f"Gene: {r['GENE']} | {r['HGVSp']}", ln=True)
            pdf.set_font("Helvetica", '', 10); pdf.multi_cell(0, 6, f"VAF: {r['VAF']:.2%} | DP: {r['DP']}\n")
            pdf.ln(2)
        return bytes(pdf.output())