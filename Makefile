.PHONY: clean docs

clean:
	rm src/*.cov
	rm test/*.pickle.gz

docs:
	julia --color=yes --project docs/make.jl