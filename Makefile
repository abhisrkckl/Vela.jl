.PHONY: clean docs

clean:
	rm src/*.cov
	rm test/*.pickle.gz

docs:
	julia --color=yes --project=docs docs/make.jl