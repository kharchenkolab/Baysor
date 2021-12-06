JULIA_ARCHIVE=julia-1.6.4-linux-x86_64.tar.gz
JULIA_DIR=$(shell pwd)/julia-1.6.4/
JULIA_DEPOT_DIR=$(shell pwd)/.julia/
JULIA_DEPOT_CONFIG=$(JULIA_DIR)/etc/julia/startup.jl
JULIA_EXEC_FLAG=$(JULIA_DIR)/bin/julia_success
JULIA_EXEC=$(JULIA_DIR)/bin/julia
JULIA_DEPS_FLAG = $(JULIA_DEPOT_DIR)/deps_success
JULIA_SYS_IMAGE=$(JULIA_DIR)/../BaysorSysimage.so
SEGMENTATION_EXECUTABLE=./baysor

all: $(SEGMENTATION_EXECUTABLE)
clean:
	rm -rf $(JULIA_ARCHIVE) $(JULIA_DIR) $(JULIA_DEPOT_DIR) $(SEGMENTATION_EXECUTABLE)

$(SEGMENTATION_EXECUTABLE): $(JULIA_SYS_IMAGE)
	$(JULIA_EXEC) --sysimage=$(JULIA_DIR)/../BaysorSysimage.so -e ""
	printf '#! $(JULIA_EXEC) --sysimage=$(JULIA_DIR)/../BaysorSysimage.so\n\nENV["JULIA_DEPOT_PATH"]="$(JULIA_DEPOT_DIR)"\nimport Baysor\nBaysor.run_cli()' >> $(SEGMENTATION_EXECUTABLE)
	chmod +x $(SEGMENTATION_EXECUTABLE)

$(JULIA_SYS_IMAGE): $(JULIA_DEPS_FLAG)
	$(JULIA_EXEC) -e 'ENV["JULIA_DEPOT_PATH"]="$(JULIA_DEPOT_DIR)"; using PackageCompiler; import Pkg; Pkg.activate(".."); Pkg.instantiate(); import Baysor; create_sysimage(:Baysor; precompile_execution_file="$$(dirname(pathof(Baysor)))/../bin/precompile.jl", sysimage_path="BaysorSysimage.so")'

$(JULIA_DEPS_FLAG): $(JULIA_DEPOT_CONFIG) $(JULIA_EXEC_FLAG)
	$(JULIA_EXEC) -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git")); Pkg.add(Pkg.PackageSpec(name="PackageCompiler", version="2.0.2"));'
	$(JULIA_EXEC) -e 'using Pkg; Pkg.build(); using Baysor' && touch $(JULIA_DEPS_FLAG)

$(JULIA_DEPOT_CONFIG): $(JULIA_EXEC_FLAG)
	mkdir -p $(JULIA_DEPOT_DIR)
	echo 'empty!(DEPOT_PATH); push!(DEPOT_PATH, "$(JULIA_DEPOT_DIR)")' >> $(JULIA_DEPOT_CONFIG)

$(JULIA_EXEC_FLAG): $(JULIA_ARCHIVE)
	rm -rf $(JULIA_DIR)
	tar xf $(JULIA_ARCHIVE) && touch $(JULIA_EXEC_FLAG)

$(JULIA_ARCHIVE):
	wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.4-linux-x86_64.tar.gz
