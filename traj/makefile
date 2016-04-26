CC=g++
FLAGS=-O3

all: traj_msd traj_gr traj_Fs traj_water traj_lwotp

traj_msd:
	mkdir -p ./bin
	$(CC) $(FLAGS) ./src/traj_msd.cc -o ./bin/traj_msd

traj_gr:
	mkdir -p ./bin
	$(CC) $(FLAGS) ./src/traj_gr.cc -o ./bin/traj_gr

traj_Fs:
	mkdir -p ./bin
	$(CC) $(FLAGS) ./src/traj_Fs.cc -o ./bin/traj_Fs

traj_water:
	mkdir -p ./bin
	$(CC) $(FLAGS) ./src/traj_water.cc -o ./bin/traj_water

traj_lwotp:
	mkdir -p ./bin
	$(CC) $(FLAGS) ./src/traj_lwotp.cc -o ./bin/traj_lwotp

clean: 
	rm ./bin/traj*
	rmdir bin
