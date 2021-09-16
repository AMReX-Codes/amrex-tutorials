add_executable(HeatEquation_EX0 main.cpp)

add_library(MyFunc myfunc.H)

target_link_libraries(HeatEquation_EX0 PRIVATE MyFunc)
