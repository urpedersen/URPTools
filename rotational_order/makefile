all: 
	mkdir -p bin
	g++ -O3 -Wall ./src/Cell_list.h ./src/Cell_list.cc ./src/split.h ./src/split.cc ./src/Rotational_order.h ./src/Rotational_order.cc ./src/Q6.cc -lboost_iostreams -o ./bin/Q6

test:
	time ./bin/Q6 -i examples/liq_fcc.xyz -Q0.3,1.0

profile:
	g++ -pg -O1 ./src/Cell_list.h ./src/Cell_list.cc ./src/Rotational_order.h ./src/Rotational_order.cc ./src/Q6.cc -lboost_iostreams -o ./bin/Q6_gprof
	./bin/Q6_gprof -i examples/crystallite.atom
	gprof ./bin/Q6_gprof gmon.out > gprof.txt
	# The below line assume that gprof2dot i installed
	gprof bin/Q6_gprof | ../../gprof2dot/gprof2dot.py | dot -Tpng -o profile.png

clean: 
	rm -f Q6.xyz node_connections.dat largest_cluster.xyz 
	rm -f gmon.out gprof.txt  profile.png 
	rm -f ./bin/Q6_gprof
	rm -f ./bin/Q6
	rmdir ./bin
