import os, re, glob
from typing import List, Dict, Tuple
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt
import seaborn as sns

class VCFProcessor:
    """Classe de extração e filtragem genômica"""

    def get_csq_fields(self, vcf_path: str) -> List[str]:
        """
        Objetivo: Mapear subcampos da anotação CSQ.
        Entrada: vcf_path (str).
        Lógica: Regex no header buscando 'Format:'.
        Saída: List[str] com os campos do VEP.
        """
        with open(vcf_path, 'r') as f:
            for l in f:
                if 'ID=CSQ' in l and 'Format:' in l:
                    return re.search(r'Format: (.*)\"', l).group(1).split('|')
        return []

    def get_vaf(self, met: Dict) -> float:
        """
        Objetivo: Calcular frequência alélica.
        Entrada: met (Dict) da amostra.
        Lógica: Prioriza AF; senão, razão AD (Alt/Total).
        Saída: float (0.0 a 1.0).
        """
        if 'AF' in met: return float(met['AF'].split(',')[0])
        ad = [float(x) for x in met.get('AD', '0,0').split(',')]
        return ad[1] / sum(ad) if sum(ad) > 0 else 0.0

    def get_sub_type(self, r: str, a: str) -> str:
        """
        Objetivo: Identificar tipo de substituição.
        Entrada: r (Ref), a (Alt).
        Lógica: Mapeia para par de pirimidina padrão.
        Saída: str (ex: 'C>T').
        """
        m = {'G>A': 'C>T', 'G>C': 'C>G', 'G>T': 'C>A', 'A>C': 'T>G', 'A>G': 'T>C', 'A>T': 'T>A'}
        return m.get(f"{r.upper()}>{a.upper()}", f"{r.upper()}>{a.upper()}")

    def validate(self, ann: Dict, dp: int, vaf: float, p: Dict) -> bool:
        """
        Objetivo: Filtragem multicritério.
        Entrada: ann, dp, vaf, p (Parâmetros).
        Lógica: Check Gene, Qualidade e Patogenicidade.
        Saída: bool.
        """
        g_ok, q_ok = ann.get('SYMBOL') in p['genes'], (dp >= p['dp_min'] or vaf >= p['vaf_min'])
        c_sig = ann.get('CLIN_SIG', '').lower()
        c_ok = any(x in c_sig for x in ['pathogenic', 'risk_factor']) if p['only_p'] else True
        return g_ok and q_ok and c_ok

    def parse_line(self, line: str, f: List, p: Dict, sid: str) -> List[Dict]:
        """
        Objetivo: Processar linha individual.
        Entrada: line, f (fields), p (params), sid.
        Lógica: Extrai métricas e valida contra CSQ.
        Saída: List[Dict] de variantes filtradas.
        """
        c = line.strip().split('\t')
        if c[6] != 'PASS' or 'CSQ=' not in c[7]: return []
        mt = dict(zip(c[8].split(':'), c[9].split(':')))
        dp, vaf = int(mt.get('DP', 0)), self.get_vaf(mt)
        for t in c[7].split('CSQ=')[1].split(';')[0].split(','):
            an = dict(zip(f, t.split('|')))
            if self.validate(an, dp, vaf, p):
                pos = an.get('Protein_position', '').split('/')[0]
                return [{'SAMPLEID': sid, 'GENE': an.get('SYMBOL'), 'VAF': vaf, 'DP': dp, 'TYPE': an.get('Consequence', '').split('&')[0], 
                         'SUB': self.get_sub_type(c[3], c[4]), 'POS': int(pos) if pos.isdigit() else 0, 'HGVSp': an.get('HGVSp', 'N/A'), 'CLIN': an.get('CLIN_SIG', 'N/A')}]
        return []

    def process_file(self, task: Tuple[str, Dict]) -> List[Dict]:
        """
        Objetivo: Processar um arquivo VCF.
        Entrada: task (path, p).
        Lógica: Itera parse_line em arquivo aberto.
        Saída: List[Dict] da amostra.
        """
        path, p = task
        fields, res, sid = self.get_csq_fields(path), [], os.path.basename(path).split('.')[0]
        if not fields: return []
        with open(path, 'r', encoding='utf-8') as f:
            for l in f:
                if not l.startswith('#'): res.extend(self.parse_line(l, fields, p, sid))
        return res

    def run_parallel(self, paths: List[str], p: Dict) -> List[Dict]:
        """
        Objetivo: Execução em Multiprocessing.
        Entrada: paths, p.
        Lógica: Pool de processos para acelerar leitura.
        Saída: Lista consolidada de variantes.
        """
        with ProcessPoolExecutor() as executor:
            tasks = [(f, p) for f in paths]
            raw = list(executor.map(self.process_file, tasks))
        return [item for sublist in raw for item in sublist]