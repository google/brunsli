OS := $(shell uname)
BINDIR = bin
OBJDIR = $(BINDIR)/obj
BROTLI_DIR = third_party/brotli
BROTLI_INCLUDE = $(BROTLI_DIR)/c/include
BRUNSLI_INCLUDE = c/include
SOURCES = $(wildcard c/common/*.cc) $(wildcard c/dec/*.cc) \
          $(wildcard c/enc/*.cc)
TOOLS_SOURCES = $(wildcard c/tools/*.cc)
OBJECTS = $(addprefix $(OBJDIR)/, $(SOURCES:.cc=.o))
ALL_OBJECTS = $(OBJECTS) $(addprefix $(OBJDIR)/, $(TOOLS_SOURCES:.cc=.o))
CBRUNSLI = $(addprefix $(OBJDIR)/, c/tools/cbrunsli.o)
DBRUNSLI = $(addprefix $(OBJDIR)/, c/tools/dbrunsli.o)
DIRS = $(OBJDIR)/c/common $(OBJDIR)/c/dec $(OBJDIR)/c/enc \
       $(OBJDIR)/c/tools
LIBBROTLI=brotli
CFLAGS += -O2 -std=c++11 -ffunction-sections
LDFLAGS += -Wl,-gc-sections
ifeq ($(os), Darwin)
  CPPFLAGS += -DOS_MACOSX
endif

ifneq ($(strip $(CROSS_COMPILE)), )
	CXX=$(CROSS_COMPILE)-gcc++
	ARCH=$(firstword $(subst -, ,$(CROSS_COMPILE)))
	BROTLI_WRAPPER="qemu-$(ARCH) -L /usr/$(CROSS_COMPILE)"
endif

# The arm-linux-gnueabi compiler defaults to Armv5. Since we only support Armv7
# and beyond, we need to select Armv7 explicitly with march.
ifeq ($(ARCH), arm)
	CFLAGS += -march=armv7-a -mfloat-abi=hard -mfpu=neon
endif

all: cbrunsli dbrunsli
	@:

.PHONY: all clean libbrotli cbrunsli dbrunsli

$(DIRS):
	mkdir -p $@

$(LIBBROTLI):
	$(MAKE) -C $(BROTLI_DIR) lib

cbrunsli: $(OBJECTS) $(CBRUNSLI) $(LIBBROTLI)
	$(CXX) $(LDFLAGS) $(OBJECTS) $(CBRUNSLI) \
        -lm -L$(BROTLI_DIR) -l$(LIBBROTLI) -o $(BINDIR)/cbrunsli

dbrunsli: $(OBJECTS) $(DBRUNSLI) $(LIBBROTLI)
	$(CXX) $(LDFLAGS) $(OBJECTS) $(DBRUNSLI) \
        -lm -L$(BROTLI_DIR) -l$(LIBBROTLI) -o $(BINDIR)/dbrunsli

clean:
	rm -rf $(BINDIR)
	$(MAKE) -C third_party/brotli clean

.SECONDEXPANSION:
$(ALL_OBJECTS): $$(patsubst %.o,%.cc,$$(patsubst $$(OBJDIR)/%,%,$$@)) | $(DIRS)
	$(CXX) $(CFLAGS) $(CPPFLAGS) -I$(BROTLI_INCLUDE) -I$(BRUNSLI_INCLUDE) \
        -c $(patsubst %.o,%.cc,$(patsubst $(OBJDIR)/%,%,$@)) -o $@
