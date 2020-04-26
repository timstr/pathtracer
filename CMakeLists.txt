cmake_minimum_required(VERSION 3.10.0)

project(pathtracer CXX)


set(SFML_STATIC_LIBRARIES TRUE)
find_package(SFML 2.4 COMPONENTS system window graphics audio main REQUIRED)

if (NOT SFML_FOUND)
	message(FATAL_ERROR "Could not find SFML. Please install SFML.")
endif()

set(pathtracer_headers
    include/PathTracer.hpp
)

set(pathtracer_srcs
    src/PathTracer.cpp
)

add_executable(pathtracer PathTracerDemo.cpp ${pathtracer_headers} ${pathtracer_srcs})

target_link_libraries(pathtracer
	PUBLIC sfml-system
	PUBLIC sfml-window
	PUBLIC sfml-graphics
	PUBLIC sfml-audio
	PUBLIC sfml-main
)

set_property(TARGET pathtracer PROPERTY CXX_STANDARD 17)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT pathtracer)


if(MSVC)
	# increase warning level
	add_definitions(/W4)

	# turn warning C4715 ("not all control paths return a value") into an error
	add_definitions(/we4715)

	# disable warning C4250 ("class X inherits virtual member function Y via dominance")
	add_definitions(/wd4250)

	# turn warning C4239 (non-standard extension that allows temporaries to be bound to
	# non-const references, yay microsoft) into an error
	# This extension is without question the source of many bugs and disastrous surprises
	add_definitions(/we4239)

    # Conformance mode (to heck with all silly Microsoft language extensions)
    add_definitions(/permissive-)
else()
	SET(CMAKE_CXX_FLAGS -pthread)
endif()