cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

include (TestBigEndian)

TEST_BIG_ENDIAN(BIGENDIAN)

if (BIGENDIAN)
    message(STATUS "*** Compiling in big endian architecture - WORDS_BIGENDIAN=1")
    add_definitions( -DWORDS_BIGENDIAN=1 )
else ()
    message(STATUS "*** Compiling in little endian architecture")
	#add_definitions( -DWORDS_BIGENDIAN=0 )
endif ()

add_library(StringEncoders modp_b64.c modp_b64.h modp_b64_data.h)
