set(CMAKE_CXX_STANDARD 17)
find_package(Threads REQUIRED)

link_libraries(${GMP_LIB})
add_definitions(${GMP_DEFINIONS})


if(USE_ASM)
    if(ARCH MATCHES "arm64")
            add_definitions(-DUSE_ASM -DARCH_ARM64)
    elseif(ARCH MATCHES "x86_64")
            add_definitions(-DUSE_ASM -DARCH_X86_64)
    endif()
endif()

if(USE_ASM AND ARCH MATCHES "x86_64")

    if (CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
        set(NASM_FLAGS -fmacho64 --prefix _)
    else()
        set(NASM_FLAGS -felf64 -DPIC)
    endif()

    add_custom_command(OUTPUT ${CMAKE_SOURCE_DIR}/build/fq_asm.o
        COMMAND nasm ${NASM_FLAGS} fq.asm -o fq_asm.o
        DEPENDS ${CMAKE_SOURCE_DIR}/build/fq.asm
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

    add_custom_command(OUTPUT ${CMAKE_SOURCE_DIR}/build/fr_asm.o
        COMMAND nasm ${NASM_FLAGS} fr.asm -o fr_asm.o
        DEPENDS ${CMAKE_SOURCE_DIR}/build/fr.asm
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build)
endif()

set(FR_SOURCES
    ../build/fr.hpp
    ../build/fr.cpp
    )

if(USE_ASM)
    if(ARCH MATCHES "arm64")

       set(FR_SOURCES ${FR_SOURCES} ../build/fr_raw_arm64.s ../build/fr_raw_generic.cpp ../build/fr_generic.cpp)
    elseif(ARCH MATCHES "x86_64")
        set(FR_SOURCES ${FR_SOURCES} ../build/fr_asm.o)
    endif()
else()
    set(FR_SOURCES ${FR_SOURCES} ../build/fr_generic.cpp ../build/fr_raw_generic.cpp)
endif()

add_library(fr STATIC ${FR_SOURCES})
set_target_properties(fr PROPERTIES POSITION_INDEPENDENT_CODE ON)

link_libraries(fr)

set(FQ_SOURCES
    ../build/fq.hpp
    ../build/fq.cpp
    )

if(USE_ASM)
    if(ARCH MATCHES "arm64")
        set(FQ_SOURCES ${FQ_SOURCES} ../build/fq_raw_arm64.s ../build/fq_raw_generic.cpp ../build/fq_generic.cpp)
    elseif(ARCH MATCHES "x86_64")
        set(FQ_SOURCES ${FQ_SOURCES} ../build/fq_asm.o)
    endif()
else()
    set(FQ_SOURCES ${FQ_SOURCES} ../build/fq_raw_generic.cpp ../build/fq_generic.cpp)
endif()

add_library(fq STATIC ${FQ_SOURCES})
set_target_properties(fq PROPERTIES POSITION_INDEPENDENT_CODE ON)

link_libraries(fq)

set(LIB_SOURCES
    alt_bn128.hpp
    alt_bn128.cpp
    binfile_utils.hpp
    binfile_utils.cpp
    curve.hpp
    exp.hpp
    f2field.hpp
    fft.hpp
    fileloader.hpp
    fullprover.hpp
    fullprover.cpp
    groth16.hpp
    groth16.cpp
    logger.hpp
    logger.cpp
    misc.hpp
    misc.cpp
    multiexp.hpp
    naf.hpp
    random_generator.hpp
    scope_guard.hpp
    spinlock.hpp
    splitparstr.hpp
    wtns_utils.hpp
    zkey_utils.hpp
    )

if(USE_LOGGER)
    set(LIB_SOURCES ${LIB_SOURCES} logger.cpp)
    add_definitions(-DUSE_LOGGER)
endif()

get_filename_component(INSTALL_PREFIX_ABS "${CMAKE_INSTALL_PREFIX}" ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}")

include_directories(
    ../src
    ../build
    ../depends/json/single_include
    ../package/include/
    ${INSTALL_PREFIX_ABS}/include
    )

link_directories(
    ../package/lib/
    ${INSTALL_PREFIX_ABS}/lib
    )


add_library(rapidsnarkStatic STATIC ${LIB_SOURCES})
set_target_properties(rapidsnarkStatic PROPERTIES OUTPUT_NAME rapidsnark)

add_library(rapidsnarkStaticFrFq STATIC ${LIB_SOURCES} ${FQ_SOURCES} ${FR_SOURCES})
add_dependencies(rapidsnarkStaticFrFq fr fq)
set_target_properties(rapidsnarkStaticFrFq PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(rapidsnarkStaticFrFq PROPERTIES OUTPUT_NAME rapidsnark-fr-fq)

target_link_libraries(rapidsnarkStatic)

add_library(rapidsnark SHARED ${LIB_SOURCES})

target_link_libraries(rapidsnark tbb Threads::Threads)


if(OpenMP_CXX_FOUND)
    if(TARGET_PLATFORM MATCHES "android")
        target_link_libraries(rapidsnark -static-openmp -fopenmp)
    endif()
endif()
