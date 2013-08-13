OPTS = -I ../deps/eigen
all: common do
	do help

common: common.cc
	g++ ${OPTS} common.cc -o common -D_COMMON_TEST 

do: common.cc
	g++ ${OPTS} common.cc -o do
