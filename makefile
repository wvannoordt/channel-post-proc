run: main
	./program

main:
	mpicxx -O3 -std=c++20 -I${CMF}/include -I${PTL}/include main.cc -o program -L${CMF}/lib -lcmf -L${PTL}/lib -lPTL

clean:
	rm -f program
