all: ../output/fig%.png

../input:
	mkdir $@

../output:
	mkdir $@

../input/Project.toml ../input/Manifest.toml: ../input/%.toml: ../../setup_environment/output/%.toml | ../input
	ln -s $< $@

../output/fig%.png: Assignment1.jl DFS1977functions.jl
	julia $<
