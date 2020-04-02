JULIA_ARCHIVE=julia-1.3.1-linux-x86_64.tar.gz
JULIA_DIR=$(shell pwd)/julia-1.3.1/
JULIA_DEPOT_DIR=$(shell pwd)/.julia/
JULIA_DEPOT_CONFIG=$(JULIA_DIR)/etc/julia/startup.jl
JULIA_EXEC_FLAG=$(JULIA_DIR)/bin/julia_success
JULIA_EXEC=$(JULIA_DIR)/bin/julia
JULIA_DEPS_FLAG = $(JULIA_DEPOT_DIR)/deps_success
JULIA_SYS_IMAGE=$(JULIA_DIR)/lib/julia/sys.so
SEGMENTATION_EXECUTABLE=./baysor

all: $(SEGMENTATION_EXECUTABLE)
clean:
	rm -rf $(JULIA_ARCHIVE) $(JULIA_DIR) $(JULIA_DEPOT_DIR) $(SEGMENTATION_EXECUTABLE)

$(SEGMENTATION_EXECUTABLE): $(JULIA_SYS_IMAGE)
	echo '#! $(JULIA_EXEC)\n\nENV["JULIA_DEPOT_PATH"]="$(JULIA_DEPOT_DIR)"\nimport Baysor\nBaysor.run_cli()' >> $(SEGMENTATION_EXECUTABLE)
	chmod +x $(SEGMENTATION_EXECUTABLE)

$(JULIA_SYS_IMAGE): $(JULIA_DEPS_FLAG)
	$(JULIA_EXEC) -e 'using GR; using PackageCompiler; import Baysor; ENV["JULIA_DEPOT_PATH"]="$(JULIA_DEPOT_DIR)"; bp=dirname(pathof(Baysor)) * "/../"; ENV["LD_LIBRARY_PATH"]=dirname(pathof(GR)) * "/../deps/gr/lib/"; create_sysimage(:Baysor; project=bp, replace_default=true)'

$(JULIA_DEPS_FLAG): $(JULIA_DEPOT_CONFIG) $(JULIA_EXEC_FLAG)
	$(JULIA_EXEC) -e 'using Pkg; Pkg.add(PackageSpec(url="git@github.com:hms-dbmi/Baysor.git")); Pkg.add(["Plots", "PackageCompiler", "GR"]);'
	$(JULIA_EXEC) -e 'using Pkg; Pkg.build(); using Baysor' && touch $(JULIA_DEPS_FLAG)

$(JULIA_DEPOT_CONFIG): $(JULIA_EXEC_FLAG)
	mkdir -p $(JULIA_DEPOT_DIR)
	echo 'empty!(DEPOT_PATH); push!(DEPOT_PATH, "$(JULIA_DEPOT_DIR)")' >> $(JULIA_DEPOT_CONFIG)

$(JULIA_EXEC_FLAG): $(JULIA_ARCHIVE)
	rm -rf $(JULIA_DIR)
	tar xf $(JULIA_ARCHIVE) && touch $(JULIA_EXEC_FLAG)

$(JULIA_ARCHIVE):
	wget https://julialang-s3.julialang.org/bin/linux/x64/1.3/julia-1.3.1-linux-x86_64.tar.gz
