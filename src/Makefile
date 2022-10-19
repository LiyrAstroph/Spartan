all:disk spec obs

disk:const.f specfun.f diskvar.f diskrad.f disksol.f diskset.f datagen.f diskmain.f diskflux.f
	gfortran -o disk diskmain.f
spec:const.f specfun.f specvar.f speccomp.f speccal.f specset.f
	gfortran -o spec specset.f
obs:obsspec.f obsvar.f raytrac.f
	gfortran -o obs obsspec.f

clean:
	rm -rf disk spec obs	
