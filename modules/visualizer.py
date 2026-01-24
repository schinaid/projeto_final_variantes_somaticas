import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import List

class BioVisualizer:
    """Gerenciador de visualizações genômicas para dashboards e laudos técnicos."""

    def plot_gene_bar(self, df: pd.DataFrame):
        """
        Objetivo: Gerar gráfico de barras de frequência mutacional.
        Entrada: df (pd.DataFrame) -> Variantes filtradas.
        Lógica: Countplot ordenado pela contagem decrescente de genes.
        Saída: Objeto Figure do Matplotlib.
        """
        fig, ax = plt.subplots(figsize=(10, 5))
        sns.countplot(data=df, x='GENE', palette='viridis', 
                      order=df['GENE'].value_counts().index, ax=ax)
        plt.title("Frequência Mutacional por Gene", fontweight='bold')
        plt.xticks(rotation=45)
        return fig

    def plot_risk_pie(self, df: pd.DataFrame, files: List[str]):
        """
        Objetivo: Proporção de risco da coorte (Amostras mutadas vs Total).
        Entrada: df (variantes), files (lista de nomes de arquivos VCF).
        Lógica: Compara IDs únicos com variantes vs total de arquivos na pasta.
        Saída: Objeto Figure do Matplotlib.
        """
        mutated = set(df['SAMPLEID'].unique())
        sids = [f.split('.')[0] for f in files]
        st_series = pd.Series(['SIM' if s in mutated else 'NÃO' for s in sids])
        fig, ax = plt.subplots()
        st_series.value_counts().plot.pie(autopct='%1.1f%%', 
                                          colors=['#ff7675', '#55efc4'], ax=ax)
        plt.title("Proporção de Risco", fontweight='bold')
        return fig

    def plot_oncoprint(self, df: pd.DataFrame):
        """
        Objetivo: Matriz OncoPrint (Paisagem Mutacional).
        Entrada: df (pd.DataFrame).
        Lógica: Pivot table binária mapeando presença de variantes por amostra.
        Saída: Objeto Figure do Matplotlib.
        """
        mtx = df.pivot_table(index='GENE', columns='SAMPLEID', 
                             values='TYPE', aggfunc='first').fillna('WT')
        fig, ax = plt.subplots(figsize=(12, len(mtx)*0.4 + 2))
        sns.heatmap(mtx.applymap(lambda x: x != 'WT'), 
                    cmap=['#f2f2f2', '#6c5ce7'], cbar=False, linewidths=1, ax=ax)
        plt.title("OncoPrint: Paisagem Mutacional", fontweight='bold')
        return fig

    def plot_lollipop(self, df: pd.DataFrame, gene: str):
        """
        Objetivo: Visualizar hotspots de mutação na proteína.
        Entrada: df (pd.DataFrame), gene (str).
        Lógica: Stem plot cruzando posição proteica e VAF (Variant Allele Frequency).
        Saída: Objeto Figure do Matplotlib.
        """
        sub = df[df['GENE'] == gene]
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.stem(sub['POS'], sub['VAF'], basefmt=" ")
        sns.scatterplot(data=sub, x='POS', y='VAF', hue='TYPE', s=100, ax=ax)
        plt.title(f"Lollipop Plot: Distribuição em {gene}", fontweight='bold')
        return fig