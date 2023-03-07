set( DEC_NAME thrplay )
set( DEC2_NAME h264dec )

option( ENABLE_DISPLAY "Enables Display using SDL2" ON )

find_package( SDL2 )

if( SDL2_FOUND )
	add_executable( ${DEC_NAME}
		threadplayer.cpp unithread.h md5.c md5.h my_getopt.c my_getopt.h getopt.h frames.h filewrite.h m2decoder.h
	)
	target_include_directories( ${DEC_NAME} PUBLIC ${PROJECT_BINARY_DIR}/src/app . ../lib ${SDL2_INCLUDE_DIRS} )
	target_link_libraries( ${DEC_NAME} libm2dec ${SDL2_LIBRARIES} )
else( SDL2_FOUND )
	set( ENABLE_DISPLAY OFF )
endif( SDL2_FOUND )

add_executable( ${DEC2_NAME}
	h264dec.cpp lowlevel.c md5.c md5.h my_getopt.c my_getopt.h getopt.h frames.h filewrite.h m2decoder.h
)
target_include_directories( ${DEC2_NAME} PUBLIC ${PROJECT_BINARY_DIR}/src/app . ../lib ${SDL2_INCLUDE_DIRS} )
target_link_libraries( ${DEC2_NAME} libm2dec ${SDL2_LIBRARIES} )

configure_file(config.h.cmake config.h)
