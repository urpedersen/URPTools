all: src/Q6.cc
	mkdir -p bin
	g++ -O3 ./src/Cell_list.h ./src/Cell_list.cc ./src/Rotational_order.h ./src/Rotational_order.cc ./src/Q6.cc -lboost_iostreams -o ./bin/Q6

clean: 
	rm ./bin/Q6
	rmdir ./bin
