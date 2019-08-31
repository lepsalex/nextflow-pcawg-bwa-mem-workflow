FROM ubuntu:18.04

# Install.
RUN \
    apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y wget

WORKDIR /usr/local/bin

RUN \
    wget https://dl.min.io/client/mc/release/linux-amd64/mc \
    && chmod +x mc

# Define default command.
CMD ["bash"]