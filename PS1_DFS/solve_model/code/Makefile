all: ../output/DFS1977_example_a.txt ../output/DFS1977_example_b.txt

../input:
	mkdir $@

../input/Project.toml ../input/Manifest.toml: ../input/%.toml: ../../setup_environment/output/%.toml | ../input
	ln -s $< $@

../output/DFS1977_example_%.txt: Assignment1.jl
	julia $<