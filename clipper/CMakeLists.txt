set(SRCS clipper.hpp clipper.cpp)
__add_real_target(clipper lib SOURCE ${SRCS})

INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include/ FILES_MATCHING PATTERN "*.hpp")