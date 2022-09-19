# To build: docker build -t hihisiv-webapp .

FROM debian
RUN apt-get update && apt-get -y install python3 python3-pip python3-psycopg2 python3-pandas && apt-get clean
RUN pip3 install --no-cache-dir streamlit pyvis
RUN useradd -ms /bin/bash streamlit
RUN sed -i 's/MESSAGE_SIZE_LIMIT = 50 \* int(1e6)/MESSAGE_SIZE_LIMIT = 50 \* int(1e7)/g' /usr/local/lib/python3.9/dist-packages/streamlit/server/server_util.py
USER streamlit
WORKDIR /home/streamlit
COPY --chown=streamlit:streamlit HIHISIV_logo.png hihisiv_webapp.py /home/streamlit/
COPY --chown=streamlit:streamlit .streamlit /home/streamlit/.streamlit
EXPOSE 8501
CMD streamlit run hihisiv_webapp.py
