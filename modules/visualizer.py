import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import List

class BioVisualizer:
    """
    Objetivo: Transformar dados tabulares em evidências visuais interpretáveis.
    Lógica: Utiliza Seaborn e Matplotlib para gerar gráficos de coorte e individuais.
    """

    def plot_gene_bar(self, df: pd.DataFrame):
        """
        Objetivo: Frequência mutacional por gene.
        Lógica: Gráfico de barras ordenado pela incidência de variantes.
        Saída: Figure (Matplotlib).
        """
        fig, ax = plt.subplots(figsize=(10, 5))
        sns.countplot(data=df, x='GENE', palette='viridis', 
                      order=df['GENE'].value_counts().index, ax=ax)
        plt.title("Prevalência Genômica na Coorte", fontweight='bold')
        plt.xticks(rotation=45)
        return fig

    def plot_risk_pie(self, df: pd.DataFrame, files: List[str]):
        """
        Objetivo: Proporção de risco da coorte.
        Lógica: Compara amostras positivas vs negativas.
        Saída: Figure (Matplotlib).
        """
        mutated = set(df['SAMPLEID'].unique())
        sids = [f.split('.')[0] for f in files]
        st = pd.Series(['SIM' if s in mutated else 'NÃO' for s in sids])
        fig, ax = plt.subplots()
        st.value_counts().plot.pie(autopct='%1.1f%%', colors=['#ff7675', '#55efc4'], ax=ax)
        plt.title("Status de Risco da Amostras", fontweight='bold')
        return fig

    def plot_oncoprint(self, df: pd.DataFrame):
        """
        Objetivo: OncoPrint (Matriz Mutacional).
        Lógica: Pivot table binária mapeando variantes por paciente.
        Saída: Figure (Matplotlib).
        """
        mtx = df.pivot_table(index='GENE', columns='SAMPLEID', 
                             values='TYPE', aggfunc='first').fillna('WT')
        fig, ax = plt.subplots(figsize=(12, len(mtx)*0.4 + 2))
        # Heatmap binário: Mutado (Purple) vs Wild-Type (Grey)
        sns.heatmap(mtx.applymap(lambda x: x != 'WT'), 
                    cmap=['#f2f2f2', '#6c5ce7'], cbar=False, ax=ax)
        return fig

    def plot_lollipop(self, df: pd.DataFrame, gene: str):
        """
        Objetivo: Lollipop Plot (Hotspots proteicos).
        Lógica: Distribuição espacial de mutações na proteína vs VAF.
        Saída: Figure (Matplotlib).
        """
        sub = df[df['GENE'] == gene]
        fig, ax = plt.subplots(figsize=(10, 4))
        # Eixo X = Posição Aminoácido | Eixo Y = Frequência Alélica
        ax.stem(sub['POS'], sub['VAF'], basefmt=" ")
        sns.scatterplot(data=sub, x='POS', y='VAF', hue='TYPE', s=100, ax=ax)
        return fig

    def plot_signatures(self, df: pd.DataFrame):
        """
        Objetivo: Assinatura Mutacional (Substituições).
        Lógica: Perfil de trocas químicas de bases nitrogenadas.
        Saída: Figure (Matplotlib).
        """
        fig, ax = plt.subplots(figsize=(8, 4))
        order = ['C>A','C>G','C>T','T>A','T>C','T>G']
        sns.countplot(data=df, x='SUB', palette='Set2', order=order, ax=ax)
        plt.title("Assinaturas: Perfil de Substituição de DNA", fontweight='bold')
        return fig