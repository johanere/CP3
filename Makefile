output: main.o vec3.o celestialbody.o system.o newtoniangravity.o solver.o
	c++ main.o vec3.o celestialbody.o system.o newtoniangravity.o solver.o -o output

test.o: main.cpp
	c++ -c main.cpp

vec3.o: vec3.cpp vec3.h
	c++ -c vec3.cpp

celestialbody.o: celestialbody.cpp celestialbody.h
	c++ -c celestialbody.cpp

system.o: system.cpp system.h
	c++ -c system.cpp

newtoniangravity.o: newtoniangravity.cpp newtoniangravity.h
	c++ -c newtoniangravity.cpp

solver.o: solver.cpp solver.h
	c++ -c solver.cpp


clean:
	rm *.o output
