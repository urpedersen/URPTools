# cluster_analysis


## Build
Just type
```
make
```
This program has no library dependencies.

## Usage
The program expects a file named node_connections.dat that should look something like
```
10    # Example with 10 nodes enumerated 0 to 9
0 1   # with nodes connections as follows:
1 2   #    8-7-6
1 8   #    |   |
3 5   #  0-1   3-5  4-9
8 7   #    |
7 6   #    2
6 3   # Notes: Characters after #'s are comments. The first line gives the number of nodes.
4 9   #        The following lines list connections between nodes.
```

### Usage example
```
make
cd Examples
../bin/cluster_analysis
```
