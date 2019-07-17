FROM pinellolab/stream

RUN conda install -c anaconda pyyaml -y

COPY main.py definition.yml /code/

RUN chmod +x /code/*

WORKDIR /code

ENTRYPOINT ["/code/main.py"]
