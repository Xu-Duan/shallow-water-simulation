compiler = gfortran
target = swan
f90s = helper.f90 timestep.f90 swan.f90
os = $(patsubst %.f90, %.o, $(f90s))

$(target):$(os)
	$(compiler) $(os) -o $(target)

%.o:%.f90
	$(compiler) -c $< -o $@

clean:
	rm -rf *.o *.mod *.out

cleanall:
	rm -rf *.o *.mod *.out $(target) output