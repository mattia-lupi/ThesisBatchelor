CXX = g++
# -MMD -MP: Generates dependency files (.d)
# -Iinclude: Looks for headers in the include/ folder
CXXFLAGS = -std=c++14 -O2 -I/opt/homebrew/include -Iinclude -Iplt -MMD -MP
LDFLAGS = -L/opt/homebrew/lib -larmadillo

TARGET = main.out

# Source files: all .cpp in root (main.cpp) AND all .cpp in src/
SRCS = $(wildcard *.cpp) $(wildcard src/*.cpp)
# Object files list (main.o src/funzioni.o etc.)
OBJS = $(SRCS:.cpp=.o)
# Dependency files list
DEPS = $(OBJS:.o=.d)

all: $(TARGET)

# Link everything
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Generic rule: Works for both main.cpp and src/file.cpp
# because the structure is preserved in $(OBJS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	time ./$(TARGET)

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)
