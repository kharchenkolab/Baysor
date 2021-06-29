FROM julia:latest

RUN apt-get update && apt-get install -y build-essential

## Jupyter

RUN apt-get install -y python3 python3-pip vim

RUN pip3 install jupyterlab numpy scipy matplotlib seaborn pandas sklearn scikit-image \
        jupyter_contrib_nbextensions jupyter_nbextensions_configurator jupyterthemes

RUN pip3 install -Iv six==1.12.0

RUN julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.build(); using IJulia;'

### jupyter notebook --no-browser --port=8989 --ip=0.0.0.0 --allow-root ./

## Julia Baysor envitonment

RUN julia -e 'using Pkg; Pkg.add("PackageCompiler")'
RUN julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/hms-dbmi/Baysor.git")); Pkg.build();'

RUN julia -e 'import Baysor, Pkg; Pkg.activate(dirname(dirname(pathof(Baysor)))); Pkg.instantiate();'
RUN julia -e 'using PackageCompiler; import Baysor, Pkg; Pkg.activate(dirname(dirname(pathof(Baysor)))); create_sysimage(:Baysor; precompile_execution_file="$(dirname(pathof(Baysor)))/../bin/precompile.jl", replace_default=true, cpu_target="generic")'

RUN \
    printf "#!/usr/bin/env julia\n\nimport Baysor\nBaysor.run_cli()" >> /bin/baysor && \
    chmod +x /bin/baysor

ENTRYPOINT ["/bin/bash"]
WORKDIR /root/
