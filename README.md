# trabalho_final_variantes_somaticas
# Pipeline de Classificação de Risco Genômico em Mielofibrose (MF)
```
Este projeto consiste em uma pipeline de bioinformática avançada para a classificação, filtragem e visualização de variantes somáticas de alto risco em pacientes com Mielofibrose.
```
🖥️ Interface do Usuário (Dashboard)

O dashboard foi projetado para fornecer uma visão clara da maturidade dos dados e da paisagem mutacional da coorte em tempo real.

Principais Seções:
🏠 Visão Geral: Exibe a frequência de mutações por gene, a proporção de risco da coorte e uma tabela consolidada de amostras.

🧩 OncoPrint: Mapa de calor binário que ilustra a distribuição de variantes em todos os pacientes simultaneamente.

🍭 Lollipop Plot: Visualização espacial das mutações na estrutura da proteína para identificar hotspots.

📊 Assinaturas: Perfil de substituição de bases nitrogenadas para análise de processos mutagênicos.

📄 PDF: Gerador de laudos técnicos automáticos com evidências gráficas e metadados de filtragem.

O pipeline é totalmente conteinerizado com **Docker**, permitindo que a análise seja reprodutível em qualquer ambiente sem a necessidade de instalar dependências de bioinformática localmente.

## 📋 Objetivo da Análise
Classificar cada amostra como `MAIOR_RISCO = SIM/NÃO` baseando-se na presença de pelo menos uma variante que cumpra:
1.  **Gene:** Presente no painel de alto risco (Ex: TP53, EZH2, etc).
2.  **Qualidade:** Filtro `PASS` no VCF.
3.  **Efeito Funcional:** Impacto `MODERATE/HIGH` ou consequências específicas (missense, stop_gained, etc).
4.  **Profundidade/Frequência:** DP ≥ 20 ou VAF ≥ 5%.
```
---

## 🏗️ Estrutura do Projeto

```text
analise_vcf/
├── inputs/               # Jogue seus arquivos VCF (e subpastas) aqui
├── outputs/              # Onde os arquivos .tsv serão salvos
│   └── plots/            # Subpasta para os gráficos gerados
├── .env                  # Configurações de genes e thresholds
├── docker-compose.yml    # Orquestração do container
├── Dockerfile            # Receita da imagem do sistema
├── main.py               # Código-fonte principal (Altamente comentado)
├── requirements.txt      # Dependências de bibliotecas
└── README.md             # Documentação de uso
```

⚙️ Configuração (.env)
Você pode ajustar os critérios de filtragem diretamente no arquivo .env sem alterar o código:

```
Variável	Descrição	Exemplo
GENES_ALTO_RISCO	Lista de genes alvo	TP53,EZH2,IDH1...
DP_MIN	Profundidade mínima de leitura	20
VAF_MIN	Frequência Alélica mínima	0.05
```

🚀 Como Executar
1. Pré-requisitos
Docker instalado.

Docker Compose instalado.

2. Preparação
Coloque seus arquivos VCF anotados (hg38 + VEP) dentro da pasta /inputs.

3. Execução
No terminal, dentro da pasta do projeto, execute:

Variável,Descrição,Exemplo/Padrão
GENES_ALTO_RISCO,Lista de genes alvo para o painel,"TP53,EZH2,CBL,U2AF1..."
DP_MIN,Profundidade mínima de leitura,20
VAF_MIN,Frequência Alélica (VAF) mínima,0.05
IMPACTOS_INTERESSE,Níveis de impacto VEP considerados,"MODERATE,HIGH"


TRABALHO BASEADO EM <a href="https://ashpublications.org/bloodadvances/article/5/5/1442/475395/Genomic-analysis-of-primary-and-secondary"> Genomic analysis of primary and secondary myelofibrosis redefines the prognostic impact of ASXL1 mutations: a FIM study</a>