FROM ubuntu

RUN apt-get update && apt-get -y upgrade && \
	apt-get install -y build-essential wget \
		libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev bowtie2 perl default-jdk seqtk unzip python3 python3-pip && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# samtools install
WORKDIR /tmp
ENV SAMTOOLS_VERSION=1.16.1
RUN wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 && \
	tar jxf samtools-$SAMTOOLS_VERSION.tar.bz2 && \
	rm samtools-$SAMTOOLS_VERSION.tar.bz2 && \
	cd samtools-$SAMTOOLS_VERSION && \
	./configure --prefix /usr/src/samtools-$SAMTOOLS_VERSION && \
	make && make install && \
    rm -rf /tmp/samtools-$SAMTOOLS_VERSION

ENV PATH=${PATH}:/usr/src/samtools-$SAMTOOLS_VERSION/bin

# bcftools install
WORKDIR /tmp
ENV BCFTOOLS_VERSION=1.16
RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
    tar --bzip2 -xf bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
    rm bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
    cd bcftools-$BCFTOOLS_VERSION && \
    ./configure --prefix=/usr/src/bcftools-$BCFTOOLS_VERSION && \
    make && make install && \
    rm -rf /tmp/bcftools-$BCFTOOLS_VERSION

ENV PATH=${PATH}:/usr/src/bcftools-$BCFTOOLS_VERSION/bin

WORKDIR /app

# bbmap install
RUN wget -c "https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download" -O - | tar -xz -C /app
ENV PATH=${PATH}:/app/bbmap

# trimmomatic install
RUN wget "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip" -O Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    rm Trimmomatic-0.39.zip

COPY ./requirements.txt /app
RUN pip3 install -r /app/requirements.txt

COPY . /app

ENTRYPOINT [ "/bin/bash" ]