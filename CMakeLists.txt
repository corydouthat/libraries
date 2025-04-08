cmake_minimum_required(VERSION 3.8)

project(libraries LANGUAGES CXX)
#set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

if (NOT TARGET libraries)
	add_library(libraries INTERFACE)
	#	array_list.hpp
	#	array_mat.hpp
	#	file_help.hpp
	#	linked_list.hpp
	#	mat.hpp
	#	math_misc.hpp
	#	quat.hpp
	#	vec.hpp
	#)
endif()

target_include_directories(libraries INTERFACE ${CMAKE_CURRENT_LIST_DIR})
