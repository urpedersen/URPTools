CC=g++
FLAGS=-O3 -Wall

all: traj_msd traj_gr traj_water traj_lwotp

traj_msd:
	mkdir -p ./bin
	$(CC) $(CFLAGS) ./src/traj_msd.cpp -o ./bin/traj_msd

traj_gr:
	mkdir -p ./bin
	$(CC) $(CFLAGS) ./src/traj_gr.cpp -o ./bin/traj_gr

traj_water:
	mkdir -p ./bin
	$(CC) $(CFLAGS) ./src/traj_water.cpp -o ./bin/traj_water

traj_lwotp:
	mkdir -p ./bin
	$(CC) $(CFLAGS) ./src/traj_lwotp.cpp -o ./bin/traj_lwotp

clean: 
	rm ./bin/traj*
	rmdir bin
