# Projeto Final – Classificação de Variantes Somáticas (MF)

Projeto de conclusão de curso – **Hospital Israelita Albert Einstein**.  

Este projeto implementa um pipeline para **classificação de risco prognóstico em Mielofibrose (MF)** com base na presença de **variantes somáticas em genes de alto risco**, a partir de arquivos VCF anotados com VEP (hg38).

A implementação segue as regras propostas na atividade do curso e é baseada no artigo:

> *Genomic analysis of primary and secondary myelofibrosis*  
> Blood Advances, 2021  
> https://ashpublications.org/bloodadvances/article/5/5/1442/475395/Genomic-analysis-of-primary-and-secondary

---

## Objetivo

Classificar cada amostra como:

- **MAIOR_RISCO = SIM / NÃO**

com base na presença de variantes somáticas patogênicas em genes de alto risco associados a pior prognóstico em MF.

Além disso, gerar uma flag independente:

- **TP53_PRESENTE = SIM / NÃO**

---

## Entrada de Dados

- **30 arquivos VCF**
- Genoma de referência: **hg38**
- Um VCF por amostra
- Arquivos já:
  - liftOver (hg19 → hg38)
  - anotados com **VEP**
- Diretório de entrada: `inputs/`

### Genes avaliados (painel fixo)

- **TP53**
- **EZH2**
- **CBL**
- **U2AF1**
- **SRSF2**
- **IDH1**
- **IDH2**
- **NRAS**
- **KRAS**

---

## Regra de Classificação

Uma amostra é considerada **MAIOR_RISCO = SIM** se existir **≥ 1 variante** que cumpra **todas** as condições abaixo:

### Critérios da variante

- `GENE` (ou `SYMBOL`) pertence ao painel definido
- `FILTER = PASS`
- Variante com **efeito funcional**, atendendo a pelo menos um critério:
  - `IMPACT ∈ {MODERATE, HIGH}` **ou**
  - `Consequence` contém:
    - `missense_variant`
    - `stop_gained`
    - `frameshift_variant`
    - `splice_`
    - `start_lost`
- Critério técnico:
  - `DP ≥ 20` **ou**
  - `VAF ≥ 0.05` (5%)

### Flag adicional

- **TP53_PRESENTE = SIM / NÃO**
  - Mesma regra acima, porém filtrando exclusivamente variantes no gene **TP53**

---

## Fluxo do Pipeline

1. Leitura dos 30 arquivos VCF
2. Extração do cabeçalho VEP (campo `CSQ`)
3. Parsing das variantes aprovadas (`FILTER = PASS`)
4. Cálculo de métricas:
   - DP
   - VAF
   - Tipo de substituição
5. Aplicação das regras biológicas e técnicas
6. Consolidação dos resultados por variante e por amostra
7. Geração de tabelas finais e visualizações
8. Exportação de gráficos e PDFs técnicos

---

## Estrutura do Projeto
```
projeto_final_variantes_somaticas/
├── app.py 
├── Dockerfile 
├── docker-compose.yml 
├── requirements.txt 
├── .env 
├── spawn.sh 
├── spawn.bat 
│ 
├── inputs/ 
│ ├── *.vcf 
│ └── leiam-me.docx 
│ 
├── modules/ 
│ ├── processor.py 
│ ├── visualizer.py 
│ └── reporter.py 
│ 
└── outputs/ 
├── variants_high_risk.tsv 
├── sample_risk.tsv 
└── plots/ 
└── dashboard_resultados.png 
```
---

## Descrição dos Diretórios

### `inputs/`

Contém os **VCFs de entrada**, um por amostra, já anotados com VEP e prontos para processamento.

---

### `modules/`

Código modularizado do pipeline.

#### `processor.py`
Responsável pelo **processamento genômico**:

- Leitura paralela dos VCFs
- Parsing do campo `CSQ`
- Cálculo de VAF
- Classificação de substituições
- Aplicação das regras de filtragem biológica
- Geração da lista de variantes de alto risco

#### `visualizer.py`
Responsável pelas **visualizações**:

- Frequência de mutações por gene
- Status de risco da coorte
- OncoPrint (matriz mutacional)
- Lollipop plots (hotspots proteicos)
- Assinaturas mutacionais

#### `reporter.py`
Responsável pela **geração de relatórios PDF**:

- Inserção de parâmetros de filtragem
- Inclusão de gráficos
- Listagem detalhada das variantes por amostra

---

### `outputs/`

Resultados finais do pipeline.

#### `variants_high_risk.tsv`
Uma linha por variante filtrada, contendo no mínimo:

- `SAMPLEID`
- `CHROM`
- `POS`
- `REF`
- `ALT`
- `GENE`
- `Consequence`
- `IMPACT`
- `FILTER`
- `DP`
- `VAF`

#### `sample_risk.tsv`
Uma linha por amostra, contendo:

- `SAMPLEID`
- `MAIOR_RISCO` (SIM/NÃO)
- `TP53_PRESENTE` (SIM/NÃO)
- `GENES_ALTO_RISCO_ENCONTRADOS`
- `N_VARIANTES_ALTO_RISCO`

#### `plots/`
Gráficos gerados automaticamente pelo dashboard.

---

## Aplicação (`app.py`)

Interface interativa construída com **Streamlit**, permitindo:

- Definição dinâmica de thresholds (DP, VAF, gnomAD)
- Execução do pipeline
- Visualização dos resultados por abas:
  - Visão geral da coorte
  - OncoPrint
  - Lollipop
  - Assinaturas mutacionais
  - Exportação de PDF por amostra

---

## Docker e Execução

### Dockerfile

- Base: `python:3.9-slim`
- Instala dependências do sistema
- Instala bibliotecas Python
- Executa o Streamlit na porta `8501`

### docker-compose.yml

- Cria o serviço `vcf-analyzer`
- Mapeia volumes:
  - `inputs/` → `/app/inputs`
  - `outputs/` → `/app/outputs`
- Usa variáveis definidas no `.env`

### `.env`

Define parâmetros do ambiente, incluindo:

- Diretório de entrada
- Diretório de saída
- Lista de genes de alto risco

---

### Execução

Linux / macOS:
```bash
./spawn.sh
```

Windows:
```
spawn.bat
```

Ou diretamente:
```
docker compose up --build
```

Resumo

Este projeto implementa um pipeline reprodutível, modular e containerizado para análise de variantes somáticas em MF, atendendo integralmente aos critérios técnicos e científicos exigidos na atividade final do curso do Einstein.




