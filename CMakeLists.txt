cmake_minimum_required(VERSION 2.8.3)
project(solver_comparison)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)
# CMake module(s) path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_STANDARD 11)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS roscpp rospy roslib)
find_package(IPOPT REQUIRED)
find_package(SNOPT REQUIRED)
FIND_PACKAGE(ACADO REQUIRED)
find_package(Ceres REQUIRED)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES mpc_acado
#  CATKIN_DEPENDS roscpp rospy
#  DEPENDS system_lib
)

include_directories(
# include
  ${catkin_INCLUDE_DIRS}
)
INCLUDE_DIRECTORIES(. ${ACADO_INCLUDE_DIRS})

add_executable(acado_node src/acado.cpp)
target_link_libraries(acado_node ${catkin_LIBRARIES} ${ACADO_SHARED_LIBRARIES} ${LINKER_LIBS})

add_executable(ipopt_node src/ipopt.cpp)
target_link_libraries(ipopt_node ${catkin_LIBRARIES} ${IPOPT_LIBRARIES} ${LINKER_LIBS})

add_executable(snopt_node src/snopt.cpp)
target_link_libraries(snopt_node ${catkin_LIBRARIES} ${SNOPT_LIBRARIES} ${LINKER_LIBS})

add_executable(ros_test_node src/ros_test.cpp)
target_link_libraries(ros_test_node ${catkin_LIBRARIES} ${LINKER_LIBS})