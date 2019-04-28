CXXFLAGS := -std=c++17 -pedantic -pedantic-errors -Wall -Wextra -Werror
LDFLAGS := -lstdc++fs

OBJ_DIR := objs
SRC := $(wildcard *.cpp)

all: build raiden

OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

$(OBJ_DIR)/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: all build clean debug release

build:
	@mkdir -p $(OBJ_DIR)

raiden: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) -o $@ 

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O3
release: all

clean:
	@rm -rvf $(OBJ_DIR)/*
	@rm -rvf raiden
