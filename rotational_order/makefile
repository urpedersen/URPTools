all: src/Q6.cc
	mkdir -p bin
	g++ -O3 -Wall ./src/Cell_list.h ./src/Cell_list.cc ./src/Rotational_order.h ./src/Rotational_order.cc ./src/Q6.cc -lboost_iostreams -o ./bin/Q6

gprof:
	g++ -pg -O3 ./src/Cell_list.h ./src/Cell_list.cc ./src/Rotational_order.h ./src/Rotational_order.cc ./src/Q6.cc -lboost_iostreams -o ./bin/Q6_gprof
	./bin/Q6_gprof -i examples/liq_fcc.xyz -o /dev/null.xyz
	gprof ./bin/Q6_gprof gmon.out > gprof.txt

clean: 
	rm ./bin/Q6
	rmdir ./bin
