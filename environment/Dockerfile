# hash:sha256:51a878fce4b1eaea0cc37af067dc23e6827754ff0b1c55d1dce907d04c0d2724
FROM registry.apps.codeocean.com/codeocean/mambaforge3:22.11.1-4-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip3 install -U --no-cache-dir \
    pandas==2.0.3 \
    pillow \
    plotly==5.16.1 \
    python-dotenv==1.0.0 \
    scipy==1.11.2 \
    statsmodels==0.14.0 \
    streamlit==1.26.0
