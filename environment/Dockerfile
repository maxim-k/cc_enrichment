# hash:sha256:89fbdaa604a5ce90c9c8ca8773fb902d51581e47ab146d4a03ccb04078ade6fe
FROM registry.apps.codeocean.com/codeocean/mambaforge3:22.11.1-4-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip3 install -U --no-cache-dir \
    pandas==2.0.3 \
    pillow==9.5.0 \
    plotly==5.16.1 \
    python-dotenv==1.0.0 \
    scipy==1.11.2 \
    statsmodels==0.14.0 \
    streamlit==1.27.0

ADD "https://github.com/coder/code-server/releases/download/v4.9.0/code-server-4.9.0-linux-amd64.tar.gz" /.code-server/code-server.tar.gz
	
RUN cd /.code-server \
	&& tar -xvf code-server.tar.gz \
	&& rm code-server.tar.gz \
	&& ln -s /.code-server/code-server-4.9.0-linux-amd64/bin/code-server  /usr/bin/code-server
