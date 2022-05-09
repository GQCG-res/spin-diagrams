FROM gqcg/gqcp:latest

COPY . .
RUN conda install -c conda-forge plotly plotly_express python-kaleido --freeze-installed

ENTRYPOINT bash
