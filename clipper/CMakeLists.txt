set(SRCS clipper.hpp clipper.cpp)
__add_real_target(clipper lib SOURCE ${SRCS})

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(LIB_DEST lib/debug)
	set(BIN_DEST bin/debug)
else()
	set(LIB_DEST lib/release)
	set(BIN_DEST bin/release)
endif()

install(TARGETS clipper
		EXPORT clipper
		LIBRARY DESTINATION ${LIB_DEST}
		ARCHIVE DESTINATION ${LIB_DEST}
		RUNTIME DESTINATION ${BIN_DEST})

install(FILES clipper.hpp
		DESTINATION include)