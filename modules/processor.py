import os, re, glob # Importação de bibliotecas para manipulação de arquivos e expressões regulares
from typing import List, Dict, Tuple # Importação de tipos para tipagem estática
from concurrent.futures import ProcessPoolExecutor # Importação para processamento paralelo em múltiplos núcleos

class VCFProcessor:
    '''
    Descrição: Classe principal para o processamento e filtragem de variantes em Mielofibrose (MF).
    Lógica: Executa parsing paralelo de VCFs, extração de anotações CSQ e validação de risco biológico.
    '''

    def __init__(self):
        '''
        Descrição: Inicializa o painel de genes e termos de consequência biológica.
        Parâmetros: Nenhum.
        Entrada: Nenhuma.
        Saída: Instância da classe configurada.
        Lógica: Define listas estáticas baseadas no consenso médico para Mielofibrose.
        '''
        self.target_genes = ['TP53', 'EZH2', 'CBL', 'U2AF1', 'SRSF2', 'IDH1', 'IDH2', 'NRAS', 'KRAS'] # Painel de genes MF
        self.target_cons = ['missense_variant', 'stop_gained', 'frameshift_variant', 'splice_', 'start_lost'] # Consequências alvo

    def get_csq_fields(self, vcf_path: str) -> List[str]:
        '''
        Descrição: Mapeia os campos da anotação CSQ (VEP) definidos no cabeçalho do arquivo VCF.
        Parâmetros: 
            - vcf_path (str): Caminho absoluto ou relativo do arquivo VCF.
        Entrada: String contendo o caminho do arquivo.
        Saída: List[str] contendo os nomes das colunas da anotação (ex: SYMBOL, IMPACT).
        Lógica: Varre o header procurando por 'ID=CSQ' e extrai o formato via Regex para garantir o mapeamento correto.
        '''
        with open(vcf_path, 'r') as f: # Abre o arquivo VCF em modo leitura
            for line in f: # Itera sobre cada linha do arquivo
                if 'ID=CSQ' in line and 'Format:' in line: # Identifica a linha de metadado do CSQ
                    match = re.search(r'Format: (.*)\"', line) # Busca o padrão do formato entre aspas
                    return match.group(1).split('|') if match else [] # Divide os campos pelo caractere pipe
        return [] # Retorna lista vazia caso não encontre a anotação CSQ no header

    def calc_vaf(self, metrics: Dict) -> float:
        '''
        Descrição: Calcula a frequência alélica da variante ($$VAF$$).
        Parâmetros:
            - metrics (Dict): Dicionário contendo os dados da coluna SAMPLE mapeados.
        Entrada: Dicionário de metadados genômicos da linha.
        Saída: float representando a frequência entre 0.0 e 1.0.
        Lógica: Utiliza a profundidade alélica (AD) para realizar o cálculo matemático Alt / Total.
        '''
        if 'AF' in metrics: return float(metrics['AF'].split(',')[0]) # Retorna AF se fornecida diretamente
        ad_values = [float(x) for x in metrics.get('AD', '0,0').split(',')] # Converte AD para lista de floats
        ref_count, alt_count = ad_values[0], ad_values[1] # Separa contagem de Referência e Alternativa
        total_depth = ref_count + alt_count # Calcula a profundidade total na posição
        return alt_count / total_depth if total_depth > 0 else 0.0 # Aplica a fórmula matemática da VAF

    def get_sub_type(self, ref: str, alt: str) -> str:
        '''
        Descrição: Normaliza as substituições de base única (SNV) para análise de assinaturas.
        Parâmetros:
            - ref (str): Base nitrogenada de referência.
            - alt (str): Base nitrogenada alternativa.
        Entrada: Duas strings representando as bases genômicas.
        Saída: str formatada representando a troca normalizada (ex: 'C>T').
        Lógica: Mapeia as substituições para os 6 pares de pirimidina padrão para evitar redundância.
        '''
        sub_map = {'G>A': 'C>T', 'G>C': 'C>G', 'G>T': 'C>A', 'A>C': 'T>G', 'A>G': 'T>C', 'A>T': 'T>A'} # Mapa químico
        pair = f"{ref.upper()}>{alt.upper()}" # Concatena o par original em caixa alta
        return sub_map.get(pair, pair) # Retorna o mapeamento ou o próprio par se já for padrão

    def validate_variant(self, ann: Dict, dp: int, vaf: float, p: Dict) -> bool:
        '''
        Descrição: Valida se a variante atende aos critérios biológicos e técnicos de risco.
        Parâmetros:
            - ann (Dict): Dicionário com anotação CSQ processada.
            - dp (int): Profundidade de leitura (DP).
            - vaf (float): Frequência alélica calculada.
            - p (Dict): Parâmetros de thresholds (DP_min, VAF_min, gnomAD).
        Entrada: Dados da variante e dicionário de limites.
        Saída: bool indicando aprovação ou rejeição da variante.
        Lógica: Aplica lógica booleana AND/OR conforme as regras de negócio de Mielofibrose.
        '''
        pop_af = float(ann.get('gnomAD_AF') or 0) # Captura frequência populacional (0 se N/A)
        gene_ok = ann.get('SYMBOL') in self.target_genes # Verifica presença no painel MF
        impact_ok = ann.get('IMPACT') in ['HIGH', 'MODERATE'] # Verifica severidade do impacto
        cons_ok = any(c in ann.get('Consequence', '') for c in self.target_cons) # Verifica termos funcionais
        metric_ok = (dp >= p['dp_min'] or vaf >= p['vaf_min']) # Aplica regra de qualidade técnica dinâmica
        return gene_ok and (impact_ok or cons_ok) and metric_ok and pop_af <= p['max_pop_af'] # Filtro final

    def parse_line(self, line: str, fields: List, p: Dict, sid: str) -> List[Dict]:
        '''
        Descrição: Realiza o parsing de uma linha individual do VCF e extrai variantes qualificadas.
        Parâmetros:
            - line (str): Linha bruta do arquivo de dados.
            - fields (List): Lista de nomes dos campos CSQ.
            - p (Dict): Parâmetros de thresholds.
            - sid (str): Identificador da amostra.
        Entrada: String da linha, metadados e parâmetros técnicos.
        Saída: List[Dict] contendo a variante processada ou lista vazia se filtrada.
        Lógica: Decompõe a linha por tabulação, extrai métricas do SAMPLE e itera sobre transcritos CSQ.
        '''
        cols = line.strip().split('\t') # Divide a linha bruta em colunas por tabulação
        if cols[6] != 'PASS' or 'CSQ=' not in cols[7]: return [] # Ignora variantes sem PASS ou sem CSQ
        format_keys = cols[8].split(':') # Extrai as chaves do campo FORMAT
        sample_vals = cols[9].split(':') # Extrai os valores do campo SAMPLE
        metrics = dict(zip(format_keys, sample_vals)) # Mapeia metadados da amostra para dicionário
        dp, vaf = int(metrics.get('DP', 0)), self.calc_vaf(metrics) # Obtém profundidade e calcula VAF
        raw_csq = cols[7].split('CSQ=')[1].split(';')[0] # Isola a string de anotação CSQ
        for transcript in raw_csq.split(','): # Itera sobre múltiplos transcritos anotados pelo VEP
            ann = dict(zip(fields, transcript.split('|'))) # Mapeia campos do VEP para valores
            if self.validate_variant(ann, dp, vaf, p): # Submete variante ao motor de validação
                p_pos = ann.get('Protein_position', '').split('/')[0] # Extrai posição da proteína
                return [{ 'SAMPLEID': sid, 'CHROM': cols[0], 'POS': cols[1], 'REF': cols[3], 'ALT': cols[4],
                          'GENE': ann.get('SYMBOL'), 'VAF': vaf, 'DP': dp, 'TYPE': ann.get('Consequence', '').split('&')[0],
                          'SUB': self.get_sub_type(cols[3], cols[4]), 'PROT_POS': int(p_pos) if p_pos.isdigit() else 0,
                          'HGVSp': ann.get('HGVSp', 'N/A'), 'CLIN': ann.get('CLIN_SIG', 'N/A'), 'IMPACT': ann.get('IMPACT') }]
        return [] # Retorna lista vazia caso nenhum transcrito passe na validação

    def run_parallel(self, paths: List[str], p: Dict) -> List[Dict]:
        '''
        Descrição: Orquestra o processamento paralelo da coorte de arquivos VCF.
        Parâmetros:
            - paths (List): Lista de caminhos físicos dos arquivos.
            - p (Dict): Parâmetros de filtragem.
        Entrada: Lista de strings e dicionário de thresholds.
        Saída: List[Dict] consolidada de toda a análise de coorte.
        Lógica: Distribui a carga de processamento entre todos os cores lógicos via ProcessPoolExecutor.
        '''
        with ProcessPoolExecutor() as executor: # Inicializa o pool de processos paralelo
            tasks = [(f, p) for f in paths] # Cria a lista de tuplas de tarefas
            results = list(executor.map(self.process_file_worker, tasks)) # Mapeia tarefas nos cores
        return [item for sublist in results for item in sublist] # Achata a lista de listas em vetor único

    def process_file_worker(self, task: Tuple[str, Dict]) -> List[Dict]:
        '''Descrição: Worker individual para leitura. Lógica: Itera sobre linhas ignorando header.'''
        path, p = task # Desempacota o caminho e os parâmetros
        fields, final_res, sid = self.get_csq_fields(path), [], os.path.basename(path).split('.')[0] # Setup
        if not fields: return [] # Cancela processamento se o arquivo não tiver CSQ
        with open(path, 'r', encoding='utf-8') as f: # Abre o VCF para leitura
            for line in f: # Varre cada linha de dados
                if not line.startswith('#'): # Filtra linhas que não são de metadados
                    final_res.extend(self.parse_line(line, fields, p, sid)) # Acumula variantes qualificadas
        return final_res # Retorna o conjunto de variantes da amostra específica