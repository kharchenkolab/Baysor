FROM julia:latest

RUN apt-get update && apt-get install -y build-essential

## Jupyter

RUN apt-get install -y python3 python3-pip

RUN pip3 install jupyterlab numpy scipy matplotlib seaborn pandas sklearn scikit-image \
        jupyter_contrib_nbextensions jupyter_nbextensions_configurator jupyterthemes

RUN pip3 install -Iv six==1.12.0

RUN \
    jupyter nbextensions_configurator enable && \
    jupyter contrib nbextension install

RUN \
    jupyter nbextension enable codefolding/main && \
    jupyter nbextension enable codefolding/edit && \
    jupyter nbextension enable collapsible_headings/main && \
    jupyter nbextension enable notify/notify && \
    jupyter nbextension enable ruler/main && \
    jupyter nbextension enable table_beautifier/main && \
    jupyter nbextension enable comment-uncomment/main && \
    jupyter nbextension enable hide_header/main && \
    jupyter nbextension enable livemdpreview/livemdpreview && \
    jupyter nbextension enable spellchecker/main && \
    jupyter nbextension enable execution_dependencies/execution_dependencies && \
    jupyter nbextension enable freeze/main && \
    jupyter nbextension enable tree-filter/index && \
    jupyter nbextension enable init_cell/main && \
    jupyter nbextension enable move_selected_cells/main && \
    jupyter nbextension enable toc2/main

RUN julia -e 'using Pkg; Pkg.add("IJulia"); Pkg.build(); using IJulia;'

### jupyter notebook --no-browser --port=8989 --ip=0.0.0.0 --allow-root ./

## Julia Baysor envitonment

RUN julia -e 'using Pkg; Pkg.add(["Plots", "PackageCompiler", "GR"]); Pkg.add(PackageSpec(url="https://github.com/hms-dbmi/Baysor.git")); Pkg.build(); using Baysor;'
RUN julia -e 'using GR; using PackageCompiler; import Baysor; bp=dirname(pathof(Baysor)) * "/../"; ENV["LD_LIBRARY_PATH"]=dirname(pathof(GR)) * "/../deps/gr/lib/"; compile_incremental("$bp/Project.toml", "$bp/bin/build.jl", force=true)'

RUN \
    echo "#!/usr/bin/env julia\n\nimport Baysor\nBaysor.run_cli()" >> /bin/segment_data && \
    chmod +x /bin/segment_data

ENTRYPOINT ["/bin/bash"]
WORKDIR /root/