FROM julia:latest

RUN apt-get update && apt-get install -y build-essential

## Jupyter

RUN apt-get install -y python3 python3-pip vim

RUN pip3 install jupyterlab numpy scipy matplotlib seaborn pandas sklearn scikit-image

RUN pip3 install -Iv six==1.12.0

RUN julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.build(); using IJulia;'

### jupyter notebook --no-browser --port=8989 --ip=0.0.0.0 --allow-root ./

## Julia Baysor envitonment

RUN julia -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="PackageCompiler", version="2.0.6"))'
RUN julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/hms-dbmi/Baysor.git")); Pkg.build();'

RUN julia -e 'import Baysor, Pkg; Pkg.activate(dirname(dirname(pathof(Baysor)))); Pkg.instantiate();'
RUN julia -e 'using PackageCompiler; import Baysor, Pkg; Pkg.activate(dirname(dirname(pathof(Baysor)))); \
    create_sysimage(:Baysor; precompile_execution_file="$(dirname(pathof(Baysor)))/../bin/precompile.jl", sysimage_path="/root/BaysorSysimage.so")'

RUN \
      printf "#!/usr/local/julia/bin/julia --sysimage=/root/BaysorSysimage.so\n\nimport Baysor\nBaysor.run_cli()" >> /bin/baysor && \
    chmod +x /bin/baysor

RUN baysor --help

ENTRYPOINT ["/bin/bash"]
WORKDIR /root/
