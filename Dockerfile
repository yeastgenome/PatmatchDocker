# ---------- builder: clone repo ----------
FROM ubuntu:24.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Which ref to build (branch, tag, or commit SHA)
ARG GIT_REF=fix_enzymes_do_not_cut
RUN git clone https://github.com/yeastgenome/PatmatchDocker.git /tmp/PatmatchDocker \
 && cd /tmp/PatmatchDocker \
 && git checkout "${GIT_REF}"

# ---------- runtime: apache + venv ----------
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    LANG=C.UTF-8 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      apache2 \
      libapache2-mod-wsgi-py3 \
      python3 \
      python3-pip \
      python3-virtualenv \
      ca-certificates \
      net-tools \
 && rm -rf /var/lib/apt/lists/*

# App files from builder
COPY --from=builder /tmp/PatmatchDocker/www /var/www
COPY --from=builder /tmp/PatmatchDocker/FlaskApp.conf /etc/apache2/sites-available/FlaskApp.conf

# Create tmp dir, make venv (use apt virtualenv to avoid ensurepip issues), install Python deps
RUN mkdir -p /var/www/tmp \
 && chmod 1777 /var/www/tmp \
 && virtualenv -p python3 /var/www/FlaskApp/FlaskApp/venv \
 && /var/www/FlaskApp/FlaskApp/venv/bin/pip install --no-cache-dir \
      "flask>=2,<3" flask-cors boto3 requests

# Enable Apache WSGI site
RUN a2enmod wsgi \
 && a2ensite FlaskApp \
 && a2dissite 000-default || true

# Container listens on port 80 (match your ECS target group / task def)
EXPOSE 80

# Optional: not used by Apache config (we set ServerName in the conf), but harmless
ENV SERVER_NAME=localhost

# Run Apache in foreground for ECS
CMD ["apachectl", "-D", "FOREGROUND"]
