
output: main.o vec3.o celestialbody.o system.o newtoniangravity.o solver.o
	c++ -std=c++11 main.o vec3.o celestialbody.o system.o newtoniangravity.o solver.o -o output

main.o: main.cpp
	g++ -c main.cpp

vec3.o: vec3.cpp vec3.h
	c++ -std=c++11 -c vec3.cpp

celestialbody.o: celestialbody.cpp celestialbody.h
	c++ -std=c++11 -c celestialbody.cpp

system.o: system.cpp system.h
	c++ -std=c++11 -c system.cpp

newtoniangravity.o: newtoniangravity.cpp newtoniangravity.h
	c++ -std=c++11 -c newtoniangravity.cpp

solver.o: solver.cpp solver.h
	c++ -std=c++11 -c solver.cpp


clean:
	rm *.o output
