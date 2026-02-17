FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create necessary directories
RUN install -d -m 755 /var/www/app && \
    install -d -m 1777 /var/www/tmp && \
    install -d -m 755 /var/www/conf

WORKDIR /var/www/app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY www/app/*.py ./
COPY www/conf/*.json /var/www/conf/

# Environment variables
ENV DATA_DIR=/data/patmatch/ \
    RESTRICTION_DATA_DIR=/data/restriction_mapper/ \
    TMP_DIR=/var/www/tmp/ \
    CONF_DIR=/var/www/conf/ \
    S3_BUCKET=""

# Health check
HEALTHCHECK --interval=30s --timeout=5s --retries=5 \
    CMD curl -fsS http://localhost:8000/ || exit 1

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
