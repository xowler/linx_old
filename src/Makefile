all: common do
	do cts -s test/bhp.gro -f test/bhp.xtc  -o a.gz  

common: common.cc
	g++ common.cc -o common -D_COMMON_TEST 

do: common.cc
	g++ common.cc -o do -lccxtc_linux_gnu_x86_64 -I ../deps/xtc -L ../deps/xtc
