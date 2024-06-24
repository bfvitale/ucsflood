FROM debian:12-slim AS build

RUN apt update -y
ARG DEBIAN_FRONTEND=noninteractive

# Install dependencies for building Python.

RUN apt install -y \
  libbz2-dev \
  libcgal-dev \
  libffi-dev \
  libgdbm-compat-dev \
  libgdal-dev \
  libgmp-dev \
  liblzma-dev \
  libmpfr-dev \
  libreadline-dev \
  libsqlite3-dev \
  libssl-dev \
  mercurial \
  python3-pip

WORKDIR /opt/app

# Build Python.

ADD https://www.python.org/ftp/python/3.11.5/Python-3.11.5.tgz .
# COPY Python-3.11.5.tgz .
RUN tar xzf Python-3.11.5.tgz
WORKDIR /opt/app/Python-3.11.5
RUN ./configure --prefix=/opt/python-3.11.5 --enable-optimizations --with-lto --with-readline && make -j10
RUN make install

# Install Python libraries.

COPY requirements.txt .
RUN /opt/python-3.11.5/bin/pip3 --no-input install -r requirements.txt
RUN /opt/python-3.11.5/bin/pip3 install gdal==3.7.1 -f https://girder.github.io/large_image_wheels

# Build and install cgalinterp.

WORKDIR /opt/app
COPY cgalinterp.tar.gz .
RUN tar xzf cgalinterp.tar.gz
WORKDIR /opt/app/cgalinterp
RUN /opt/python-3.11.5/bin/python3 -m build
RUN /opt/python-3.11.5/bin/pip3 install dist/cgalinterp-0.0.1-cp311-cp311-linux_x86_64.whl

# Build final image.

FROM debian:12-slim AS final

COPY --from=build /opt/python-3.11.5 /opt/python-3.11.5

# Set up locale and timezone. Install utilities.
# Done in one compound RUN do reduce container layers.
# TODO: we squash layers at end; perhaps unnecessary to use one compound RUN.

RUN apt update -y \
    && apt install -y --no-install-recommends \
       apt-transport-https ca-certificates gnupg \
       curl emacs-nox grass locales libreadline8 libssl3 net-tools openssl \
       openssh-server sudo tzdata zip \
    && ln -sf /usr/share/zoneinfo/America/New_York /etc/localtime \
    && echo "America/New_York" > /etc/timezone \
    && dpkg-reconfigure -f noninteractive tzdata \
    && echo 'en_US.UTF-8 UTF-8' > /etc/locale.gen \
    && echo 'LANG="en_US.UTF-8"' > /etc/default/locale \
    && locale-gen \
    && dpkg-reconfigure --frontend=noninteractive locales \
    && update-locale LANG=en_US.UTF-8 \
    && addgroup --gid 1000 ben \
    && adduser --gid 1000 --uid 1000 \
               --gecos "Ben Vitale" ben --disabled-password \
    && usermod -aG sudo ben \
    && chown -R ben:ben /home/ben \
    && install -d -m 1777 scratch \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && install -d -m 700 -o ben -g ben /home/ben/.ssh \
    && mkdir /var/run/sshd \
    && sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

# Install Gcloud SDK / CLI.
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
         gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | \
	 tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && apt-get update -y && apt install google-cloud-cli -y

# Setup ssh.
COPY --chown=ben:ben id_rsa_docker.pub /home/ben/.ssh/authorized_keys

# Setup GCP MTU in Docker.
COPY daemon.json /etc/docker/daemon.json

# Squash to one layer.

FROM scratch
COPY --from=final / /

ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANGUAGE en_US:en  
WORKDIR /home/ben
