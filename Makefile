CC = gcc

ifndef DEBUG
 DEBUG = 1
endif
ifndef OPT
 OPT = 0
endif
ifndef FGSTATIC
 FGSTATIC = 0
endif

ifdef windir
 delcall = del
 binfix = .exe
 pathfix = $(subst /,\,$1)
 OSID = 1
else
 ifeq ($(shell uname), Linux)
  delcall = rm -f
  binfix =
  pathfix = $1
  OSID = 0
 else
  ifneq ($(shell uname), Linux)
   $(error OS not identified)
  endif
 endif
endif

ifeq ($(DEBUG),1)
 CFLAGS = -O0 -g -Wall -Wextra -Wpedantic -lm -std=c99 -fopenmp
 #add -U__STRICT_ANSI__ when compiling with MSYS
else ifeq ($(DEBUG),2)
  CFLAGS = -O0 -g -pg -Wall -Wextra -Wpedantic -lm -std=c99 -fprofile-arcs -ftest-coverage -fopenmp
else ifeq ($(DEBUG),0)
 ifeq ($(OPT),1)
  CFLAGS = -O3 -ffast-math -lm -std=c99 -fopenmp
  BFLAGS = -Wl,--strip-all
 else
  CFLAGS = -O3 -lm -std=c99 -fopenmp
  BFLAGS = -Wl,--strip-debug
 endif
endif

ifndef BFLAGS
 BFLAGS =
endif

ifndef OWNFLAGS
 OWNFLAGS = -lm
endif

ifeq ($(OSID),1)
#static compile and linking for windows and freeglut
 ifeq ($(FGSTATIC),1)
  FGLSTAT = -D FREEGLUT_STATIC
  LINKOGL = -lfreeglut_static -lglu32 -lopengl32 -lwinmm -lgdi32
 else
  FGLSTAT =
  LINKOGL = -lglut -lglu32 -lopengl32
 endif
else
 FGLSTAT =
 LINKOGL = -lGL -lGLU -lglut
endif

LFLAGS = $(LINKOGL) $(OWNFLAGS)

ODIR = obj
SDIR = src

FILES = main.c draw.c control.c auxf.c font.c ray.c alloc.c lina.c refr_data.c tests.c shapes.c color.c rot.c intersec.c fresnel.c ray_aux.c msg.c tofile.c fromfile.c tofig.c viewer.c
_OBJECTS = $(patsubst %.c, %.o, $(FILES))
_GPF = $(patsubst %.c, %.gcno, $(FILES))

OBJECTS = $(patsubst %, $(ODIR)/%, $(_OBJECTS))
GPF = $(patsubst %, $(ODIR)/%, $(_GPF))

ifeq ($(OSID),1)
 DELOBJECTS := $(patsubst %, $(ODIR)\\%, $(_OBJECTS))
 DELGPF := $(patsubst %, $(ODIR)\\%, $(_GPF))
else
 DELOBJECTS := $(OBJECTS)
 DELGPF := $(GPF)
endif

artpic: $(OBJECTS)
	$(CC) $(CFLAGS) $(BFLAGS) -o $@ $^ $(LFLAGS)

$(ODIR)/%.o: $(SDIR)/%.c $(SDIR)/artpic.h
	$(CC) $(CFLAGS) -c $< -o $@ $(FGLSTAT)

$(OBJECTS): | $(ODIR)

$(ODIR):
	mkdir $(ODIR)

clean:
	$(delcall) artpic$(binfix) $(DELOBJECTS)
ifeq ($(DEBUG),2)
	$(delcall) $(DELGPF)
endif
