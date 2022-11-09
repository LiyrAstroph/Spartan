all:disk spec obs

SRC = ./src

disk:$(SRC)/const.f $(SRC)/specfun.f $(SRC)/diskvar.f $(SRC)/diskrad.f $(SRC)/disksol.f \
     $(SRC)/diskset.f $(SRC)/datagen.f $(SRC)/diskmain.f $(SRC)/diskflux.f
	gfortran -w -o disk $(SRC)/diskmain.f
spec:$(SRC)/const.f $(SRC)/specfun.f $(SRC)/specvar.f $(SRC)/speccomp.f $(SRC)/speccal.f \
     $(SRC)/specset.f
	gfortran -w -o spec $(SRC)/specset.f
obs:$(SRC)/obsspec.f $(SRC)/obsvar.f $(SRC)/raytrac.f
	gfortran -w -o obs $(SRC)/obsspec.f

clean:
	rm -rf disk spec obs	
