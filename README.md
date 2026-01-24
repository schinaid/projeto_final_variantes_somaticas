# 🧬 Genomic Risk Analysis – Myelofibrosis (Enterprise)

Este repositório contém uma **plataforma robusta de bioinformática** para a **classificação, filtragem e visualização de variantes somáticas de alto risco** em pacientes com **Mielofibrose**.

O sistema foi desenvolvido seguindo princípios de **Engenharia de Software Enterprise**, com **arquitetura modular orientada a objetos (POO)** e otimização para **processamento paralelo em clusters de alta performance (HPC)**.

---

## 🖥️ Interface e Visualização Geral

O **dashboard principal** fornece uma visão consolidada da maturidade dos dados e da **paisagem mutacional da coorte**, permitindo uma avaliação rápida do **risco clínico**.

### 🏠 Componentes da Visão Geral

- **Frequência Mutacional por Gene**  
  Identifica genes com maior incidência de variantes aprovadas pelos filtros técnicos.

- **Proporção de Risco**  
  Gráfico que discrimina amostras de:
  - **Maior Risco (SIM)**
  - **Risco Padrão (NÃO)**  
  A classificação é baseada na presença de mutações de **alto impacto clínico**.

- **Tabela Consolidada de Amostras**  
  Visão agregada por paciente com os seguintes metadados:

| SAMPLEID              | risco_maior | genes_encontrados              | n_variantes |
|-----------------------|-------------|--------------------------------|-------------|
| liftOver_WP216_hg38   | SIM         | U2AF1, TP53                   | 2           |
| liftOver_WP280_hg38   | NÃO         | —                              | 0           |
| liftOver_WP306_hg38   | SIM         | CBL, EZH2, IDH1, U2AF1         | 6           |

---

## 🏗️ Arquitetura Modular

O software segue padrões rigorosos de **Engenharia de Software Enterprise**, garantindo:

- Separação clara entre **lógica de negócio** e **camada de visualização**
- Alta **manutenibilidade**
- Código testável e escalável

> 🔒 Cada função respeita o limite máximo de **7 linhas de lógica funcional**, reduzindo complexidade cognitiva.

### 🔧 Componentes do Sistema

- **VCFProcessor**  
  Motor de processamento paralelo responsável por:
  - Parsing de arquivos VCF
  - Filtragem biológica
  - Execução em cluster HPC

- **BioVisualizer**  
  Encapsula a geração de visualizações científicas:
  - OncoPrint
  - Lollipop Plot
  - Assinaturas Mutacionais  
  *(Matplotlib + Seaborn)*

- **ReportManager**  
  Consolida dados e evidências gráficas em **laudos técnicos em PDF**, com foco em legibilidade clínica.

- **Streamlit App**  
  Orquestra a interface do usuário e utiliza `st.session_state` para persistência de dados entre interações.

---

## 🔍 Definição de Filtros e Parâmetros Técnicos

O rigor científico do pipeline é garantido por **múltiplas camadas de filtragem configuráveis**.

### ⚙️ Métricas de Qualidade (QC)

- **DP (Depth)**  
  Profundidade mínima de leitura para garantir confiança estatística.

- **VAF (Variant Allele Frequency)**  
  Proporção de alelos mutados, essencial para inferir clonalidade:

\[
VAF = \frac{\text{Leituras Alternativas}}{\text{Total de Leituras}}
\]

---

### 🔬 Filtros Populacionais e Clínicos

- **gnomAD Max AF**  
  Remove variantes com frequência populacional acima do limiar definido (ex: > 1%), caracterizando polimorfismos germinativos comuns.

- **ClinVar – Pathogenic Only**  
  Quando ativado, exibe apenas variantes classificadas como:
  - *Patogênicas*
  - *Provavelmente Patogênicas*

---

## 📄 Relatórios Técnicos (PDF)

O sistema gera **laudos automatizados** contendo:

- **Metadados de Filtragem**  
  Tabela de referência com parâmetros de DP, VAF e gnomAD aplicados.

- **Evidências Gráficas**  
  Inclusão de **Lollipop Plot** para validação de hotspots proteicos.

- **Detalhamento Técnico**  
  Lista completa de variantes com:
  - HGVSp
  - Classificação ClinVar
  - Profundidade de leitura

---

## 🚀 Considerações Finais

Este projeto foi concebido para **ambientes clínicos e de pesquisa**, oferecendo:
- Robustez científica
- Escalabilidade computacional
- Clareza na comunicação de risco genômico

> Ideal para pipelines de **medicina de precisão**, **oncogenômica** e **pesquisa translacional**.

---
