FROM ubuntu:20.04

RUN apt-get update && apt-get install -y build-essential
RUN apt install -y curl
RUN curl -fsSL https://install.julialang.org | sh -s -- -y

ENV PATH="/root/.juliaup/bin:${PATH}"
ENV JULIA_NUM_THREADS=30

# ### Ignore cache (https://stackoverflow.com/questions/35134713/disable-cache-for-specific-run-commands)
# # ARG CACHEBUST=1
# RUN julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git#develop"));'
# RUN julia -e 'using Pkg; Pkg.develop("Baysor");'

# ENTRYPOINT ["/bin/bash"]
# WORKDIR /root/.julia/dev/Baysor/
