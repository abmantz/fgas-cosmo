# Instructions for using the source distribution with COSMOMC

Apart from compiling the package itself (separately), the installation consists of pasting a few short blocks of code into the COSMOMC Makefile and source code. These instructions work as of the July 2020 COSMOMC version and are not guaranteed in perpetuity.

1. Go to your COSMOMC directory.
2. Unzip the download and run `make` in the `fgas-cosmo/src` directory to compile the code. You may need to customize the Makefile, though it's fairly simple.
3. Copy or link `fgas-cosmo/cosmomc/fgas.ini` to one of the COSMOMC batch directories (optional, for convenience).
4. Copy or link `fgas-cosmo/cosmomc/fgas.f90` to the COSMOMC `source/` directory (required).
5. In the COSMOMC source/DataLikelihoods.f90, add these lines in the obvious places.
```
use fgas
call FgasLikelihood_Add(DataLikelihoods, Ini)
```
6. In the COSMOMC `source/Makefile`, make these additions as indicated:
```
### at the top is fine
# fgas module
FGAS ?= ../fgas-cosmo
ifneq ($(FGAS),)
FGASO = fgas.o
FGASF = fgas.f90
endif

### after the DATAMODULES definition
ifneq ($(FGAS),)
DATAMODULES += $(OUTPUT_DIR)/$(FGASO)
endif

### can go anywhere after the OBJFILES and LINKFLAGS definitions
ifneq ($(FGAS),)
OBJFILES += $(FGAS)/src/fwrapper.o
LINKFLAGS += -lstdc++ -L$(FGAS)/src -lfgas
endif
```

You should now be able to `make` COSMOMC. If you've previously compiled, you might need to do a make clean first to force CAMB to recompile completely as well.

Like any other data set, you can use the fgas data with an `INCLUDE` or `DEFAULT` statement in your ini file that calls `fgas.ini`. This module can be used with COSMOMC's "CMB" and "astro" model parametrizations (but not "background", as this lacks the baryon density).
