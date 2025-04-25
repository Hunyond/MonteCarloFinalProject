CXX = g++
CXXFLAGS = -O3 -pthread -Wall -Wextra -std=c++20 -fdiagnostics-color=always -funroll-loops -march=native

TARGET = PhotonTransport.exe
SRCS = $(wildcard *.cpp)
DEPS = $(wildcard *.hpp)

all: $(TARGET)

$(TARGET): $(SRCS) $(DEPS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS)

clean:
	$(RM) $(TARGET)