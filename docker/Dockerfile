FROM r-base:4.3.3
WORKDIR /app
RUN apt-get update && \
  apt-get install -y libxml2-dev

# install general packages 
#RUN apt-get -y install libcurl4-openssl-dev
#RUN apt-get -y install libssl-dev
#RUN apt-get -y install python3 python3-pip
RUN set -xe \
    && apt-get update \
    && apt-get -o Dpkg::Options::="--force-overwrite" install python3-pip -y
RUN python3 -m pip install --upgrade pip --break-system-packages 
RUN pip --version
RUN apt-get install -y procps

# Install required python libraries
COPY requirements-pip.txt .
RUN pip install --no-cache-dir --upgrade pip --break-system-packages  \
  && pip install --no-cache-dir -r requirements-pip.txt --break-system-packages 

# Install R packages from binary (much quicker than using install.packages)
COPY ./requirements-rbin.txt .
RUN cat requirements-rbin.txt | xargs apt-get install -y -qq

# Install any chosen R packages that aren't available as binaries
COPY ./requirements-src.R .
RUN Rscript requirements-src.R




