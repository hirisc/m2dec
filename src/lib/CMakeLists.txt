include(TestBigEndian)
test_big_endian(WORDS_BIGENDIAN)
check_type_size("int *" SIZEOF_INT_P)
configure_file(config.h.cmake config.h)

set( DEC_LIB libm2dec )
file( GLOB INC_FILES "*.h" )
file( GLOB CXX_FILES "*.cpp" )
file( GLOB C_FILES "*.c" )
add_library( ${DEC_LIB} STATIC ${C_FILES} ${CXX_FILES} )

target_include_directories( ${DEC_LIB}
	PUBLIC . ${CMAKE_BINARY_DIR}/src/lib
)
target_link_libraries( ${DEC_LIB} )
