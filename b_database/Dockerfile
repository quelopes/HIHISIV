# To build: docker build -t hihisiv-postgres .

FROM postgres
ENV POSTGRES_PASSWORD postgres
ENV POSTGRES_DB hihisiv
COPY init.sql /docker-entrypoint-initdb.d/
