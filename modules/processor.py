import os, re, glob
from typing import List, Dict, Tuple
from concurrent.futures import ProcessPoolExecutor

class VCFProcessor:
    """
    Objetivo: Motor de processamento genômico de alta performance.
    Entrada: Arquivos VCF brutos anotados.
    Lógica: Realiza o parsing paralelo e a filtragem biológica rigorosa.
    Saída: DataFrames estruturados para análise clínica.
    """

    def get_csq_fields(self, vcf_path: str) -> List[str]:
        """
        Objetivo: Mapear a estrutura da anotação CSQ (VEP) no cabeçalho do VCF.
        Entrada: vcf_path (str) -> Caminho físico do arquivo.
        Lógica: Localiza a linha ##INFO=<ID=CSQ e extrai os campos via Regex.
        Saída: List[str] com a ordem dos campos (ex: SYMBOL, CLIN_SIG).
        """
        with open(vcf_path, 'r') as f:
            for line in f:
                # Localiza a definição do formato da anotação VEP
                if 'ID=CSQ' in line and 'Format:' in line:
                    return re.search(r'Format: (.*)\"', line).group(1).split('|')
        return []

    def calc_vaf(self, met: Dict) -> float:
        """
        Objetivo: Calcular a Variant Allele Frequency (VAF).
        Entrada: met (Dict) -> Dicionário da coluna SAMPLE.
        Lógica: Prioriza campo AF; fallback para AD (Alt / Total).
        Saída: float (Frequência entre 0.0 e 1.0).
        """
        # Verifica se o campo AF já existe (frequência alélica pronta)
        if 'AF' in met: return float(met['AF'].split(',')[0])
        # Cálculo manual baseado na profundidade alélica (AD)
        ad = [float(x) for x in met.get('AD', '0,0').split(',')]
        # Lógica: Alt / (Ref + Alt)
        return ad[1] / sum(ad) if sum(ad) > 0 else 0.0

    def get_sub_signature(self, r: str, a: str) -> str:
        """
        Objetivo: Classificar substituições para assinaturas mutacionais.
        Entrada: r (Ref), a (Alt) -> Bases nitrogenadas.
        Lógica: Normaliza substituições para o par de pirimidina (C/T).
        Saída: str formatada (ex: 'C>T').
        """
        # Mapeamento para as 6 substituições padrão de Watson-Crick
        m = {'G>A': 'C>T', 'G>C': 'C>G', 'G>T': 'C>A', 'A>C': 'T>G', 'A>G': 'T>C', 'A>T': 'T>A'}
        pair = f"{r.upper()}>{a.upper()}"
        return m.get(pair, pair)

    def validate(self, ann: Dict, dp: int, vaf: float, p: Dict) -> bool:
        """
        Objetivo: Validar se a variante atende critérios biológicos e técnicos.
        Entrada: ann (Anotação), dp, vaf, p (Thresholds da UI).
        Lógica: Check de Gene, Qualidade, gnomAD (populacional) e ClinVar.
        Saída: bool (Verdadeiro se aprovada).
        """
        gnom_af = float(ann.get('gnomAD_AF') or 0)
        # Filtros: Presença no painel, DP/VAF mínimos e frequência populacional
        g_ok, q_ok = ann.get('SYMBOL') in p['genes'], (dp >= p['dp_min'] or vaf >= p['vaf_min'])
        # Filtro de patogenicidade clínica opcional
        sig = ann.get('CLIN_SIG', '').lower()
        c_ok = any(x in sig for x in ['pathogenic', 'risk_factor']) if p['only_p'] else True
        return g_ok and q_ok and c_ok and gnom_af <= p['max_pop_af']

    def parse_line(self, line: str, fields: List, p: Dict, sid: str) -> List[Dict]:
        """
        Objetivo: Processar uma linha e extrair variantes qualificadas.
        Entrada: line, fields (CSQ), p (Parâmetros), sid (Sample ID).
        Lógica: Parsing de colunas, extração de métricas e validação biológica.
        Saída: List[Dict] contendo a variante processada.
        """
        c = line.strip().split('\t')
        # Filtra apenas variantes PASS com anotação CSQ presente
        if c[6] != 'PASS' or 'CSQ=' not in c[7]: return []
        mt = dict(zip(c[8].split(':'), c[9].split(':')))
        dp, vaf = int(mt.get('DP', 0)), self.calc_vaf(mt)
        # Itera sobre transcritos anotados
        for t in c[7].split('CSQ=')[1].split(';')[0].split(','):
            an = dict(zip(fields, t.split('|')))
            if self.validate(an, dp, vaf, p):
                pos = an.get('Protein_position', '').split('/')[0]
                return [{'SAMPLEID': sid, 'GENE': an.get('SYMBOL'), 'VAF': vaf, 'DP': dp, 
                         'SUB': self.get_sub_signature(c[3], c[4]), 'POS': int(pos) if pos.isdigit() else 0,
                         'TYPE': an.get('Consequence', '').split('&')[0], 'HGVSp': an.get('HGVSp', 'N/A'), 'CLIN': an.get('CLIN_SIG', 'N/A')}]
        return []

    def process_file(self, task: Tuple[str, Dict]) -> List[Dict]:
        """
        Objetivo: Worker paralelo para processar um arquivo VCF.
        Entrada: task (Tupla contendo caminho e parâmetros).
        Lógica: Abre arquivo, identifica campos e itera linhas de dados.
        Saída: List[Dict] consolidada da amostra.
        """
        path, p = task
        flds, res, sid = self.get_csq_fields(path), [], os.path.basename(path).split('.')[0]
        if not flds: return []
        with open(path, 'r', encoding='utf-8') as f:
            for l in f:
                # Ignora metadados (linhas com #)
                if not l.startswith('#'): res.extend(self.parse_line(l, flds, p, sid))
        return res

    def run_parallel(self, paths: List[str], p: Dict) -> List[Dict]:
        """
        Objetivo: Orquestrar execução paralela no cluster (EKS/Docker).
        Entrada: Lista de arquivos e dicionário de parâmetros.
        Lógica: Distribui processamento entre os núcleos da CPU.
        Saída: List[Dict] mestre de toda a coorte.
        """
        with ProcessPoolExecutor() as executor:
            tasks = [(f, p) for f in paths]
            raw = list(executor.map(self.process_file, tasks))
        return [item for sublist in raw for item in sublist]