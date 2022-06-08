FROM ubuntu:20.04 as builder

RUN DEBIAN_FRONTEND=noninteractive apt-get update \
    && apt-get install -y git wget \
    && apt-get autoremove

WORKDIR /tools

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz \
    && tar zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz \
    && rm ncbi-blast-2.13.0+-x64-linux.tar.gz \
    && git clone https://github.com/yeastgenome/BlastDocker.git

#####

FROM ubuntu:20.04 

WORKDIR /tools/ncbi-blast-2.13.0+
COPY --from=builder /tools/ncbi-blast-2.13.0+ .

WORKDIR /tools
RUN ln -s ncbi-blast-2.13.0+ blast \
    && DEBIAN_FRONTEND=noninteractive apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get upgrade -y \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        apache2 \
        libapache2-mod-wsgi-py3 \
        net-tools \
        python3-pip \
    && pip3 install Flask \
    && pip3 install -U flask-cors \
    && pip3 install virtualenv 

WORKDIR /var/www
COPY --from=builder /tools/BlastDocker/www .

WORKDIR /etc/apache2/sites-available
COPY --from=builder /tools/BlastDocker/FlaskApp.conf .

WORKDIR /var/www/FlaskApp/FlaskApp/venv
WORKDIR /var/www/FlaskApp/FlaskApp
RUN a2enmod wsgi \
    && a2ensite FlaskApp \
    && a2dissite 000-default \
    && virtualenv venv \
    && . venv/bin/activate

WORKDIR /

CMD ["apachectl", "-D", "FOREGROUND"]
