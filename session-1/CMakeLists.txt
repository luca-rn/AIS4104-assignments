find_package(Eigen3 CONFIG REQUIRED)

add_executable(session-1 main.cpp
        main.cpp
        ../math/main.cpp
        ../math/include/math/math.h)

target_link_libraries(session-1 PUBLIC math PRIVATE Eigen3::Eigen)
