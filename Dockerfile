FROM julia:latest

RUN apt-get update && apt-get install -y build-essential

## Jupyter

RUN apt-get install -y python3 python3-pip vim

RUN pip3 install jupyterlab numpy scipy matplotlib seaborn pandas sklearn scikit-image

RUN pip3 install -Iv six==1.12.0

RUN julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.build(); using IJulia;'

### jupyter notebook --no-browser --port=8989 --ip=0.0.0.0 --allow-root ./

## Julia Baysor envitonment

### Ignore cache (https://stackoverflow.com/questions/35134713/disable-cache-for-specific-run-commands)
ARG CACHEBUST=1
RUN julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git"));'

ENV LazyModules_lazyload false

RUN julia -e 'import Baysor, Pkg; Pkg.activate(dirname(dirname(pathof(Baysor)))); Pkg.instantiate(); Pkg.build();'
RUN echo "export PATH=/root/.julia/bin/:$PATH" >> ~/.bashrc
RUN echo "alias julia='/usr/local/julia/bin/julia --sysimage=/root/.julia/scratchspaces/cc9f9468-1fbe-11e9-0acf-e9460511877c/sysimg/libbaysor.so'" >> ~/.bashrc
RUN ln -s /root/.julia/bin/baysor /usr/local/bin/baysor

RUN /root/.julia/bin/baysor --help

ENTRYPOINT ["/bin/bash"]
WORKDIR /root/
