# 1. Base Image: Slim para eficiência de recursos
FROM python:3.9-slim

# 2. Metadados de Governança
LABEL maintainer="Anderson Alves Schinaid <andersonschinaid@hotmail.com>"
LABEL project="Análise Genômica Mielofibrose"
LABEL version="1.0.0"
LABEL description="Análise paralela de coorte de VCFs para Mielofibrose"

# 3. Variáveis de Ambiente de Sistema
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# 4. Instalação de dependências e limpeza de cache para reduzir tamanho
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    python3-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

# 5. Configuração de Diretório e Segurança
WORKDIR /app

# Criação de usuário não-root para rodar no Kubernetes (Segurança Enterprise)
RUN groupadd -r genomics && useradd -r -g genomics -m genomics \
    && chown -R genomics:genomics /app

# 6. Instalação das bibliotecas Python (Cache Layer)
COPY --chown=genomics:genomics requirements.txt .
RUN pip install -r requirements.txt

# 7. Cópia do Código Fonte com permissão correta
COPY --chown=genomics:genomics . .

# 8. Definição do usuário de execução
USER genomics

# 9. Healthcheck: Permite que o EKS monitore a saúde do Pod
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD curl --fail http://localhost:8501/_stcore/health || exit 1

# 10. Exposição da porta
EXPOSE 8501

# 11. Comando de Inicialização usando variáveis de ambiente
CMD ["streamlit", "run", "app.py"]