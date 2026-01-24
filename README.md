
# 🧬 Genomic Risk Analysis – Myelofibrosis (Enterprise)

Este repositório contém uma **plataforma de bioinformática** para a **classificação, filtragem e visualização de variantes somáticas de alto risco** em pacientes com **Mielofibrose**.

O sistema foi desenvolvido seguindo princípios de **Engenharia de Software**, com **arquitetura modular orientada a objetos (POO)** e otimização para **processamento paralelo**.

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

O software segue padrões rigorosos de **Engenharia de Software**, garantindo:

- Separação clara entre **lógica de negócio** e **camada de visualização**
- Alta **manutenibilidade**
- Código testável e escalável

> 🔒 Cada função respeita o limite máximo de **7 linhas de lógica funcional**, reduzindo complexidade cognitiva.

### 🔧 Componentes do Sistema

- **VCFProcessor**  
  Motor de processamento paralelo responsável por:
  - Parsing de arquivos VCF
  - Filtragem
  - Execução em cluster HPC

- **BioVisualizer**  
  Encapsula a geração de visualizações científicas:
  - OncoPrint
  - Lollipop Plot
  - Assinaturas Mutacionais  
  *(Matplotlib + Seaborn)*

- **ReportManager**  
  Consolida dados e evidências gráficas em **PDF**, com foco em legibilidade clínica.

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




# 🧬 Genomic Risk Analysis – Myelofibrosis (Enterprise)

Plataforma de **bioinformática clínica** para **classificação, filtragem e visualização de variantes somáticas de alto risco** em pacientes com **Mielofibrose**.

O sistema foi desenvolvido seguindo princípios rigorosos de **Engenharia de Software**, com **arquitetura modular orientada a objetos (POO)** e processamento paralelo.

---

## Como Executar e Acessar o Sistema

Após iniciar a aplicação Streamlit, acesse no navegador:

http://127.0.0.1:8501



A partir dessa interface, o usuário pode **controlar dinamicamente os filtros genômicos** e explorar diferentes visões analíticas da coorte.

---

## 🖥️ Interface Geral do Sistema

### 🔧 Painel de Controle (Sidebar)

A barra lateral funciona como o **centro de governança da análise**, permitindo ajustes em tempo real:

- **DP Min (Depth)**  
  Define a profundidade mínima de leitura para garantir confiabilidade estatística.

- **VAF Min (Variant Allele Frequency)**  
  Define a fração mínima de alelos mutados:
  
  \[
  VAF = \frac{\text{Leituras Alternativas}}{\text{Total de Leituras}}
  \]

- **gnomAD Max AF**  
  Remove variantes comuns da população saudável.  
  Valores > 0.01 (1%) indicam forte evidência de variante germinativa.

- **Apenas ClinVar Pathogenic**  
  Quando ativado, mantém apenas variantes classificadas como:
  - Patogênicas
  - Provavelmente Patogênicas

- **Iniciar Análise de VCFs**  
  Executa o pipeline completo de filtragem, agregação e visualização.

---

## 🏠 Aba 1: Visão Geral (Overview)

Esta aba funciona como o **Executive Summary da coorte**.

### O que o sistema está fazendo
Filtra variantes somáticas **em tempo real**, aplicando os thresholds de DP, VAF e gnomAD definidos na sidebar.

### O que você está visualizando

- **Gráfico de Barras – Frequência Mutacional por Gene**  
  Identifica quais genes do painel de Mielofibrose (ex: *JAK2, CALR, MPL, ASXL1*) concentram mais mutações aprovadas.

- **Gráfico de Pizza – Proporção de Risco**  
  Compara amostras classificadas como:
  - **SIM** → Presença de mutações de alto impacto
  - **NÃO** → Nenhuma variante relevante após filtragem

- **Tabela Consolidada de Amostras**  
  Lista cada `SAMPLEID`, o status de risco e os genes encontrados, facilitando auditoria e governança.

![Visão Geral](./b39298d7-01c7-4739-ad8f-3f9530a83e4b.png)

---

## 🧩 Aba 2: OncoPrint – Paisagem Mutacional

O **OncoPrint** permite visualizar a paisagem mutacional completa da coorte em uma única matriz.

### O que está fazendo
O módulo `BioVisualizer` transforma a lista de variantes em uma **matriz binária**:
- Linhas → Genes
- Colunas → Pacientes (Amostras)

### O que o gráfico explica
- Co-ocorrência de mutações
- Exclusividade mútua entre genes

### Significado
Blocos coloridos indicam presença de mutação, facilitando a identificação de **drivers genômicos centrais**.

![OncoPrint](./2b171e59-328f-4af5-bf6c-811d4061a176.png)

---

## 🍭 Aba 3: Lollipop Plot – Distribuição Proteica

Enquanto o OncoPrint analisa a coorte, o **Lollipop Plot** foca na **proteína individual**.

### O que está fazendo
Mapeia a posição genômica da variante para a coordenada proteica (`Protein_position` – VEP).

### O que você está visualizando
- **Eixo X** → Extensão da proteína
- **Eixo Y** → VAF (carga mutacional)

### Significado
Permite identificar **hotspots funcionais**, essenciais em Mielofibrose (ex: domínios quinase).

![Lollipop](./ab7032dd-0636-4de8-ad02-7218462f8da3.png)

---

## 📊 Aba 4: Assinaturas Mutacionais

Esta aba investiga **os processos biológicos subjacentes às mutações**.

### O que está fazendo
Agrupa SNVs em seis classes:
- C>A, C>G, C>T
- T>A, T>C, T>G

### O que você visualiza
Gráfico de barras com a distribuição percentual de cada substituição.

### Significado
Cada padrão reflete um processo biológico distinto.  
Por exemplo, **C>T** é fortemente associado ao envelhecimento celular.

![Assinaturas](./8d686f85-db6a-403c-b1dc-8ca7fbfe5bb8.png)

---

## 📄 Aba 5: Laudo Técnico (PDF)

Etapa final da pipeline, onde o dado vira **documento clínico ou científico oficial**.

### O que está fazendo
Integra:
- Metadados dos filtros
- Evidências gráficas (Lollipop)
- Lista detalhada de variantes

### Conteúdo do Relatório
- **Contexto de Filtro** (DP, VAF, gnomAD)
- **Interpretação Clínica (ClinVar)**
- **Detalhamento técnico (HGVSp, profundidade, impacto)**

![PDF](./1f81d343-91b0-4399-b1b5-ae8b41cdef76.png)
---

