project(sumsOfPowers)
cmake_minimum_required(VERSION 2.8)


add_definitions('-g')
add_definitions('-Wall')
add_definitions('-std=c++11')


aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})


