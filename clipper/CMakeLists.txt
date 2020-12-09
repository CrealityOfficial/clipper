set(SRCS clipper.hpp clipper.cpp clip.cpp)
__add_real_target(clipper lib SOURCE ${SRCS})
set_property(TARGET clipper PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/../)

if(INSTALL_CLIPPER)
	INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include/ FILES_MATCHING PATTERN "*.hpp")
endif()
