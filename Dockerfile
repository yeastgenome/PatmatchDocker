FROM ubuntu:24.04 as builder

RUN DEBIAN_FRONTEND=noninteractive apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y git wget \
    && DEBIAN_FRONTEND=noninteractive apt-get autoremove \
    && git clone https://github.com/yeastgenome/PatmatchDocker.git

WORKDIR /PatmatchDocker
RUN git checkout restrictionmapper_fix_6542

#####

FROM ubuntu:24.04

RUN DEBIAN_FRONTEND=noninteractive apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get upgrade -y \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        apache2 \
        libapache2-mod-wsgi-py3 \
        net-tools \
	python3-boto3 \
	python3-flask \
 	python3-flask-cors \
        python3-pip \
	python3-virtualenv \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /PatmatchDocker/www /var/www
COPY --from=builder /PatmatchDocker/FlaskApp.conf /etc/apache2/sites-available/

WORKDIR /var/www/tmp
WORKDIR /var/www/FlaskApp/FlaskApp/static
WORKDIR /var/www/FlaskApp/FlaskApp/venv
WORKDIR /var/www/FlaskApp/FlaskApp
RUN chmod 1777 /var/www/tmp \
    && a2enmod wsgi \
    && a2ensite FlaskApp \
    && a2dissite 000-default \
    && virtualenv venv \
    && . venv/bin/activate

WORKDIR /

CMD ["apachectl", "-D", "FOREGROUND"]
