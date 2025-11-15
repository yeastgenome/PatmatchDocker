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
      python3-venv \
      ca-certificates \
      net-tools \
 && rm -rf /var/lib/apt/lists/*

# Use local repo instead of git clone
COPY ./www /var/www
COPY ./FlaskApp.conf /etc/apache2/sites-available/FlaskApp.conf

RUN mkdir -p /var/www/tmp \
 && chmod 1777 /var/www/tmp \
 && python3 -m venv /var/www/FlaskApp/FlaskApp/venv \
 && /var/www/FlaskApp/FlaskApp/venv/bin/pip install --upgrade pip \
 && /var/www/FlaskApp/FlaskApp/venv/bin/pip install --no-cache-dir \
      "flask>=2,<3" flask-cors boto3 requests

RUN a2enmod wsgi \
 && a2ensite FlaskApp \
 && a2dissite 000-default || true

EXPOSE 80
ENV SERVER_NAME=localhost

CMD ["apachectl", "-D", "FOREGROUND"]
