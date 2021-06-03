JULIA_ARCHIVE=julia-1.6.0-linux-x86_64.tar.gz
JULIA_DIR=$(shell pwd)/julia-1.6.0/
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
	printf '#! $(JULIA_EXEC)\n\nENV["JULIA_DEPOT_PATH"]="$(JULIA_DEPOT_DIR)"\nimport Baysor\nBaysor.run_cli()' >> $(SEGMENTATION_EXECUTABLE)
	chmod +x $(SEGMENTATION_EXECUTABLE)

$(JULIA_SYS_IMAGE): $(JULIA_DEPS_FLAG)
	$(JULIA_EXEC) -e 'ENV["JULIA_DEPOT_PATH"]="$(JULIA_DEPOT_DIR)"; using PackageCompiler; import Baysor; create_sysimage(:Baysor; precompile_execution_file="$$(dirname(pathof(Baysor)))/../bin/precompile.jl", replace_default=true)'

$(JULIA_DEPS_FLAG): $(JULIA_DEPOT_CONFIG) $(JULIA_EXEC_FLAG)
	$(JULIA_EXEC) -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git")); Pkg.add("PackageCompiler");'
	$(JULIA_EXEC) -e 'using Pkg; Pkg.build(); using Baysor' && touch $(JULIA_DEPS_FLAG)

$(JULIA_DEPOT_CONFIG): $(JULIA_EXEC_FLAG)
	mkdir -p $(JULIA_DEPOT_DIR)
	echo 'empty!(DEPOT_PATH); push!(DEPOT_PATH, "$(JULIA_DEPOT_DIR)")' >> $(JULIA_DEPOT_CONFIG)

$(JULIA_EXEC_FLAG): $(JULIA_ARCHIVE)
	rm -rf $(JULIA_DIR)
	tar xf $(JULIA_ARCHIVE) && touch $(JULIA_EXEC_FLAG)

$(JULIA_ARCHIVE):
	wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.0-linux-x86_64.tar.gz
