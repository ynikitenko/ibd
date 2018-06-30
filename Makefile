lib : libibd.so.1

Main :
	@g++ -O3 -o $@.exe *.cpp 

libibd.so.1 :
	@g++ -Wall -fPIC -c main.cpp
	@g++ -shared -Wl,-soname,$@ -o $@ main.o
	@echo $@ done.

clean :
	@rm -f *.o

.PHONY : lib clean
# http://www.gnu.org/software/make/manual/make.html#Phony-Targets

