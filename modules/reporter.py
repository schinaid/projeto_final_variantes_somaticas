import io
from fpdf import FPDF
import pandas as pd

class ReportManager:
    """Gerenciador de laudos técnicos com foco em legibilidade e metadados clínicos."""

    def _add_metadata_section(self, pdf, p: dict):
        """
        Objetivo: Inserir tabela de referência dos parâmetros de filtragem aplicados.
        Lógica: Bloco formatado com DP min, VAF min e threshold gnomAD.
        """
        pdf.set_fill_color(240, 240, 240)
        pdf.set_font("Helvetica", 'B', 10)
        pdf.cell(0, 10, " PARÂMETROS DE FILTRAGEM (CRITÉRIOS DE CORTE)", ln=True, fill=True)
        pdf.set_font("Helvetica", '', 9)
        txt = f"DP Mín: {p['dp_min']} | VAF Mín: {p['vaf_min']:.2f} | gnomAD Max AF: {p['max_pop_af']:.3f}"
        pdf.cell(0, 8, txt, ln=True)
        pdf.ln(5)

    def _fig_to_buffer(self, fig) -> io.BytesIO:
        """Objetivo: Converter plot em buffer binário PNG para o PDF."""
        buf = io.BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
        buf.seek(0)
        return buf

    def _add_variant_block(self, pdf, r: pd.Series):
        """
        Objetivo: Listar variante com espaçamento tipográfico corrigido.
        Lógica: Usa multi_cell com altura de linha (h=6) para evitar caracteres grudados.
        """
        pdf.set_font("Helvetica", 'B', 11)
        pdf.cell(0, 8, f"Gene: {r['GENE']} | {r['HGVSp']}", ln=True)
        pdf.set_font("Helvetica", '', 10)
        info = f"VAF: {r['VAF']:.2%} | DP: {r['DP']} | ClinVar: {r['CLIN']}\n"
        pdf.multi_cell(0, 6, info)
        pdf.ln(4)

    def create_pdf(self, sid: str, df: pd.DataFrame, p: dict, figs: list = None) -> bytes:
        """
        Objetivo: Orquestrar a criação do PDF completo.
        Entrada: Amostra ID, variantes, parâmetros e lista de figuras.
        Saída: Bytes binários do documento PDF.
        """
        pdf = FPDF()
        pdf.add_page()
        pdf.set_font("Helvetica", 'B', 18)
        pdf.cell(0, 15, f"PDF GENÔMICO: {sid}", ln=True, align='C')
        self._add_metadata_section(pdf, p)
        if figs:
            for f in figs:
                pdf.image(self._fig_to_buffer(f), w=170)
                pdf.ln(5)
        [self._add_variant_block(pdf, row) for _, row in df.iterrows()]
        return bytes(pdf.output())