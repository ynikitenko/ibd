all : main lib 

lib : libibd.so.1

main : *.cpp *.h
	@g++ -O3 -o $@ *.cpp 
	@echo $@ done.

libibd.so.1 : main.cpp
	@g++ -Wall -fPIC -c main.cpp
	@g++ -shared -Wl,-soname,$@ -o $@ main.o
	@echo $@ done.

clean :
	@rm -f *.o

.PHONY : all lib clean
# http://www.gnu.org/software/make/manual/make.html#Phony-Targets

