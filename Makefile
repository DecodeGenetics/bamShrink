
CXXFLAGS+=-Iinclude
CXXFLAGS+=-pthread

# RELEASE build
CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -lz -DSEQAN_HAS_ZLIB=1 -DNDEBUG=1

# set std to c++0x to allow using 'auto' etc.
CXXFLAGS+=-std=c++0x


all: bamShrink
