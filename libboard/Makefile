CXX       = g++
DIR       = linux
CXXFLAGS  = -DUNIX_MODE -std=c++11 -pedantic -I. -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -Wno-unknown-pragmas -Wno-char-subscripts
LDFLAGS   = -O2

MESHLIB   = libboard.a
SRCS      = Board.cpp Color.cpp\
	 Path.cpp PSFonts.cpp\
	 Rect.cpp ShapeList.cpp\
	 Shapes.cpp Tools.cpp\
	 Transforms.cpp

OBJS      =  $(DIR)/Board.o $(DIR)/Color.o\
	 $(DIR)/Path.o $(DIR)/PSFonts.o\
	 $(DIR)/Rect.o $(DIR)/ShapeList.o\
	 $(DIR)/Shapes.o $(DIR)/Tools.o\
	 $(DIR)/Transforms.o
default: $(DIR)/$(MESHLIB)

$(DIR)/$(MESHLIB): $(OBJS)
	ar rcs $(DIR)/$(MESHLIB) $(OBJS) 

clean:
	rm -f $(OBJS) $(DIR)/$(MESHLIB)

$(DIR)/Board.o: Board.cpp  Board.h board/Point.h board/Rect.h board/Shapes.h board/Tools.h board/PSFonts.h
	${CXX} $(CXXFLAGS) -c Board.cpp -o $(DIR)/Board.o

$(DIR)/Color.o: Color.cpp  board/Color.h board/Tools.h
	${CXX} $(CXXFLAGS) -c Color.cpp -o $(DIR)/Color.o

$(DIR)/Path.o: Path.cpp  board/Path.h board/Transforms.h
	${CXX} $(CXXFLAGS) -c Path.cpp -o $(DIR)/Path.o

$(DIR)/PSFonts.o: PSFonts.cpp
	${CXX} $(CXXFLAGS) -c PSFonts.cpp -o $(DIR)/PSFonts.o

$(DIR)/Rect.o: Rect.cpp  board/Rect.h
	${CXX} $(CXXFLAGS) -c Rect.cpp -o $(DIR)/Rect.o

$(DIR)/ShapeList.o: ShapeList.cpp  board/ShapeList.h board/Tools.h
	${CXX} $(CXXFLAGS) -c ShapeList.cpp -o $(DIR)/ShapeList.o

$(DIR)/Shapes.o: Shapes.cpp  board/Rect.h board/Shapes.h board/Tools.h board/PSFonts.h
	${CXX} $(CXXFLAGS) -c Shapes.cpp -o $(DIR)/Shapes.o

$(DIR)/Tools.o: Tools.cpp  board/Tools.h
	${CXX} $(CXXFLAGS) -c Tools.cpp -o $(DIR)/Tools.o

$(DIR)/Transforms.o: Transforms.cpp  board/Rect.h board/Shapes.h board/ShapeList.h board/Transforms.h
	${CXX} $(CXXFLAGS) -c Transforms.cpp -o $(DIR)/Transforms.o

