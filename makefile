#maincc := main.cc
maincc := ext4fft.cc
run: main
	./extpf
	./extfft

main:
	mpicxx -O3 -std=c++20 -I${CMF}/include -I${PTL}/include main.cc    -o extpf  -L${CMF}/lib -lcmf -L${PTL}/lib -lPTL
	mpicxx -O3 -std=c++20 -I${CMF}/include -I${PTL}/include ext4fft.cc -o extfft -L${CMF}/lib -lcmf -L${PTL}/lib -lPTL

clean:
	rm -f program
