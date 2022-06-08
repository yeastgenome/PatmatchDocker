FROM ubuntu:20.04

RUN DEBIAN_FRONTEND=noninteractive apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get upgrade -y \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        apache2 \
        git \
        libapache2-mod-wsgi-py3 \
        net-tools \
        python3-pip \
    && pip3 install Flask \
    && pip3 install -U flask-cors \
    && pip3 install virtualenv \
    && pip3 install boto3 \
    && rm -rf /var/lib/apt/lists/* \
    && git clone https://github.com/yeastgenome/PatmatchDocker.git

WORKDIR /PatmatchDocker
RUN cp -pr /PatmatchDocker/www /var \
    && chmod 1777 /var/www/tmp \
    && cp -p /PatmatchDocker/FlaskApp.conf /etc/apache2/sites-available \
    && a2enmod wsgi \
    && a2ensite FlaskApp \
    && a2dissite 000-default

WORKDIR /var/www/FlaskApp/FlaskApp/static
WORKDIR /var/www/FlaskApp/FlaskApp/venv
WORKDIR /var/www/FlaskApp/FlaskApp
RUN virtualenv venv \
    && . venv/bin/activate

WORKDIR /

CMD ["apachectl", "-D", "FOREGROUND"]
