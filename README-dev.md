## Preparing a release

Update the version in the Project.toml, then:

```
export BAYSOR_VERSION=v0.6.0
# ...change the version in the Project.toml...
LazyModules_lazyload=false JULIA_NUM_THREADS=30 julia --project ./deps/build.jl app
# ...test transferability...
zip -r "baysor-x86_x64-linux-${BAYSOR_VERSION}_build.zip" LICENSE README.md ./bin/baysor/*

git push origin master
docker build -t vpetukhov/baysor:latest -t "vpetukhov/baysor:$BAYSOR_VERSION" --build-arg CACHEBUST=$(date +%s) .
git tag -a $BAYSOR_VERSION -m $BAYSOR_VERSION
git push origin master --tags

docker push -a vpetukhov/baysor
```

## Building the application with Comonicon

```julia
using Baysor
using Comonicon: Builder, Configs

options = Configs.read_options(Baysor)
Builder.build_application(Baysor, options)
```
