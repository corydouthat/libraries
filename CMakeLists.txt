cmake_minimum_required(VERSION 3.19)

project(libraries LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

# Fully guard by placing include directory inside guard
if (NOT TARGET libraries)
	add_library(libraries INTERFACE
		#including files here will force Visual Studio to show library
		src/vec.hpp
		src/quat.hpp
		src/mat.hpp
		#src/vec.cpp
		#src/quat.cpp
		#src/mat.cpp
		src/array_list.hpp
		src/array_mat.hpp
		src/linked_list.hpp
		src/math_misc.hpp
		src/file_help.hpp
	 )

	target_include_directories(libraries INTERFACE src)
endif()
