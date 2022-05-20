FROM ubuntu:20.04

RUN DEBIAN_FRONTEND=noninteractive apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get upgrade -y \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        apache2 \
        libapache2-mod-wsgi-py3 \
        net-tools \
        python3-pip \
    && pip3 install Flask \
    && pip3 install -U flask-cors \
    && pip3 install virtualenv \
    && pip3 install boto3
    && rm -rf /var/lib/apt/lists/*

WORKDIR /var/www
COPY www .
RUN chmod 1777 tmp

WORKDIR /etc/apache2/sites-available
COPY FlaskApp.conf .

WORKDIR /var/www/FlaskApp/FlaskApp/venv
WORKDIR /var/www/FlaskApp/FlaskApp

RUN a2enmod wsgi \
    && a2ensite FlaskApp \
    && a2dissite 000-default \
    && virtualenv venv \
    && . venv/bin/activate

WORKDIR /

CMD ["apachectl", "-D", "FOREGROUND"]
