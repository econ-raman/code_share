all: ../output/Project.toml ../output/Manifest.toml

../output:
	mkdir $@

../output/%.toml: setup_environment.jl
	julia $< 