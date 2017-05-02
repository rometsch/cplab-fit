RM := rm -rf

# All Target
all: Genetisch

# Tool invocations
Genetisch:
	gfortran -c ./src/pikaia.f -o ./src/pikaia.o
	g++ ./src/pikaia.o ./src/pikaia_aufruf.cpp -std=c++11 -O3 -lm -lgfortran -o Genetisch


# Other Targets
clean:
	-$(RM) ./src/pikaia.o Genetisch
	-@echo ' '
