cmake_minimum_required(VERSION 3.10.0)

project(pathtracer CXX)

set(CMAKE_CXX_STANDARD 20)

# set(SFML_STATIC_LIBRARIES TRUE)
find_package(SFML 2.4 COMPONENTS system window graphics audio REQUIRED)

if (NOT SFML_FOUND)
	message(FATAL_ERROR "Could not find SFML. Please install SFML.")
endif()

set(pathtracer_headers
    include/BasicMaterial.hpp
    include/Camera.hpp
    include/Color.hpp
    include/ColorBounce.hpp
    include/Geometry.hpp
    include/Image.hpp
    include/LinearAlgebra.hpp
    include/Object.hpp
    include/ObjectTree.hpp
    include/RandomNumberGenerator.hpp
    include/Renderer.hpp
    include/RenderSettings.hpp
    include/Scene.hpp
    include/ToneMapper.hpp
)

set(pathtracer_srcs
    src/BasicMaterial.cpp
    src/Camera.cpp
    src/Color.cpp
    src/ColorBounce.cpp
    src/Geometry.cpp
    src/Image.cpp
    src/LinearAlgebra.cpp
    src/Object.cpp
    src/ObjectTree.cpp
    src/RandomNumberGenerator.cpp
    src/Renderer.cpp
    src/RenderSettings.cpp
    src/Scene.cpp
    src/ToneMapper.cpp
)

add_library(pathtracer STATIC ${pathtracer_headers} ${pathtracer_srcs})

target_include_directories(pathtracer
	PUBLIC include
)

add_executable(pathtracerdemo PathTracerDemo.cpp)

target_link_libraries(pathtracerdemo
    PUBLIC pathtracer
	PUBLIC sfml-system
	PUBLIC sfml-window
	PUBLIC sfml-graphics
	PUBLIC sfml-audio
	#	PUBLIC sfml-main
)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT pathtracerdemo)


if(MSVC)
	# increase warning level
	add_definitions(/W4)

	# turn warning C4715 ("not all control paths return a value") into an error
	add_definitions(/we4715)

	# disable warning C4250 ("class X inherits virtual member function Y via dominance")
	add_definitions(/wd4250)

    # Sane exception handling
    # https://docs.microsoft.com/en-us/cpp/build/reference/eh-exception-handling-model?view=msvc-170
    add_definitions(/EHsc)

	# turn warning C4239 (non-standard extension that allows temporaries to be bound to
	# non-const references, yay microsoft) into an error
	# This extension is without question the source of many bugs and disastrous surprises
	add_definitions(/we4239)

    # Conformance mode (to heck with all silly Microsoft language extensions)
    add_definitions(/permissive-)

    # Speeeeeeed
    set(
        CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /fp:fast /GL"
    )

    set_property(TARGET pathtracerdemo APPEND PROPERTY LINK_FLAGS /LTCG)
else()
	SET(CMAKE_CXX_FLAGS -pthread)
endif()
