output: main.o vec3.o celestialbody.o system.o newtoniangravity.o solver.o
	g++ -std=c++11 main.o vec3.o celestialbody.o system.o newtoniangravity.o solver.o -o output

main.o: main.cpp
	g++ -std=c++11 -c main.cpp

vec3.o: vec3.cpp vec3.h
	g++ -std=c++11 -c vec3.cpp

celestialbody.o: celestialbody.cpp celestialbody.h
	g++ -std=c++11 -c celestialbody.cpp

system.o: system.cpp system.h
	g++ -std=c++11 -c system.cpp

newtoniangravity.o: newtoniangravity.cpp newtoniangravity.h
	g++ -std=c++11 -c newtoniangravity.cpp

solver.o: solver.cpp solver.h
	g++ -std=c++11 -c solver.cpp


clean:
	rm *.o output
