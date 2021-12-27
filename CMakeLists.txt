set(SRCS clipper/clipper.hpp clipper/clipper.cpp
		 clipperxyz/clipper.hpp clipperxyz/clipper.cpp)
		 
__add_real_target(clipper lib SOURCE ${SRCS}
							  INC ${CMAKE_CURRENT_SOURCE_DIR}
							  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
							  SOURCE_FOLDER
							  )

if(CC_INSTALL_CLIPPER)
	__install_directory_specif(${CMAKE_CURRENT_SOURCE_DIR}/ clipper)
endif()