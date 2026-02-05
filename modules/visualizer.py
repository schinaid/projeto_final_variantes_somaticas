import matplotlib.pyplot as plt # Importação do núcleo de plotagem Matplotlib
import seaborn as sns # Importação do Seaborn para estilização estatística
import pandas as pd # Manipulação de estruturas de dados
import os # Operações de sistema de arquivos

class BioVisualizer:
    '''
    Descrição: Gerador de evidências visuais científicas para o dashboard genômico.
    Lógica: Transforma DataFrames filtrados em representações gráficas interpretáveis.
    '''

    def plot_gene_frequency(self, df: pd.DataFrame):
        '''Descrição: Gráfico de barras de prevalência. Lógica: Countplot ordenado por incidência.'''
        if df.empty: return None # Proteção contra execução com dados nulos
        fig, ax = plt.subplots(figsize=(10, 5)) # Inicializa o container da figura
        sns.countplot(data=df, x='GENE', palette='viridis', order=df['GENE'].value_counts().index, ax=ax) # Plot
        plt.title("Prevalência Mutacional por Gene (Coorte)", fontweight='bold') # Adiciona título
        return fig # Retorna objeto Figure para renderização no App

    def plot_oncoprint(self, df: pd.DataFrame):
        '''Descrição: Matriz de Co-ocorrência. Lógica: Heatmap binário em dados pivotados.'''
        if df.empty: return None # Proteção contra execução com dados nulos
        mtx = df.pivot_table(index='GENE', columns='SAMPLEID', values='TYPE', aggfunc='first').fillna('WT') # Pivot
        fig, ax = plt.subplots(figsize=(12, len(mtx)*0.4 + 2)) # Ajusta altura dinâmica conforme nº de genes
        sns.heatmap(mtx.applymap(lambda x: x != 'WT'), cmap=['#f2f2f2', '#6c5ce7'], cbar=False, ax=ax) # Plot
        plt.title("OncoPrint: Paisagem da Coorte", fontweight='bold') # Adiciona título
        return fig # Retorna objeto Figure

    def plot_signatures(self, df: pd.DataFrame):
        '''Descrição: Perfil de Substituição de Bases. Lógica: Countplot SNV em ordem biológica.'''
        if df.empty: return None # Proteção contra execução com dados nulos
        fig, ax = plt.subplots(figsize=(8, 4)) # Inicializa container
        order = ['C>A','C>G','C>T','T>A','T>C','T>G'] # Ordem canônica de substituições
        sns.countplot(data=df, x='SUB', palette='Set2', order=order, ax=ax) # Plot
        plt.title("Assinaturas: Perfil de Substituição SNV", fontweight='bold') # Título
        return fig # Retorna objeto Figure

    def plot_risk_pie(self, df: pd.DataFrame, files: list):
        '''Descrição: Proporção de risco global. Lógica: Compara IDs mutados vs IDs totais.'''
        mutated_ids = set(df['SAMPLEID'].unique()) # IDs que possuem pelo menos uma variante de risco
        total_ids = [os.path.basename(f).split('.')[0] for f in files] # Todos os IDs processados
        risk_series = pd.Series(['SIM' if s in mutated_ids else 'NÃO' for s in total_ids]) # Mapeamento binário
        fig, ax = plt.subplots() # Inicializa container
        risk_series.value_counts().plot.pie(autopct='%1.1f%%', colors=['#ff7675', '#55efc4'], ax=ax) # Pizza
        plt.title("Status de Risco Global (Coorte)", fontweight='bold') # Título
        return fig # Retorna objeto Figure