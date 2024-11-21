CC=g++

.PHONY:
	all

libs =

FLAGS = -g -Wall -Werror

STANDARD = -std=c++14

src_dir = AAC_Project

files = $(src_dir)/AAC_Project.cpp

build_dir = bin

all: common.o graph_cycles.o AAC_Project.o
	g++ $(FLAGS) $(STANDARD) $(libs) -o $(build_dir)/AAC_Project $(build_dir)/common.o $(build_dir)/graph_cycles.o $(build_dir)/AAC_Project.o

run:
	./$(build_dir)/AAC_Project

AAC_Project.o:
	g++ $(FLAGS) $(STANDARD) $(libs) -c -o $(build_dir)/AAC_Project.o $(src_dir)/AAC_Project.cpp

common.o:
	g++ $(FLAGS) $(STANDARD) $(libs) -c -o $(build_dir)/common.o $(src_dir)/common.cpp

graph_cycles.o:
	g++ $(FLAGS) $(STANDARD) $(libs) -c -o $(build_dir)/graph_cycles.o $(src_dir)/graph_cycles.cpp