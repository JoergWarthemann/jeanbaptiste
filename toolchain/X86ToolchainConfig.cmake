set(SDK_ROOT /opt/sdks/meluxs-sdks/stm32mp1-smartboard/5.6.0)

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR ${CMAKE_HOST_SYSTEM_PROCESSOR})

set(SDK_SYSROOTS "${SDK_ROOT}/sysroots")
set(SDK_NATIVE_SYSROOT "${SDK_SYSROOTS}/x86_64-ostl_sdk-linux")

set(CMAKE_FIND_ROOT_PATH "${SDK_NATIVE_SYSROOT}")

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)
set(CMAKE_GCOV "${SDK_NATIVE_SYSROOT}/usr/bin/arm-meluxs-linux/arm-meluxs-linux-gcov")
set(CMAKE_GCOVR "${SDK_NATIVE_SYSROOT}/usr/bin/python3" "${SDK_NATIVE_SYSROOT}/usr/bin/gcovr")
set(CMAKE_XSLTPROC "${SDK_NATIVE_SYSROOT}/usr/bin/xsltproc")
set(CMAKE_XSLTSYLE "${CMAKE_CURRENT_SOURCE_DIR}/tools/gtest2html/gtest2html.xslt")

#set(CMAKE_C_COMPILER "${SDK_NATIVE_SYSROOT}/usr/bin/x86_64-ostl_sdk-linux-gcc")
#set(CMAKE_CXX_COMPILER "${SDK_NATIVE_SYSROOT}/usr/bin/x86_64-ostl_sdk-linux-g++")
#set(CMAKE_CXX_FLAGS "${ALLLANG_FLAGS} -fdiagnostics-color -ffunction-sections -fmessage-length=0 -Wall -fmacro-prefix-map=${ROOT_SOURCE_DIR}/=./" CACHE STRING "" FORCE)

set(CMAKE_EXE_LINKER_FLAGS "-pthread -lrt -ldl -lutil -ltinfo")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath=${SDK_NATIVE_SYSROOT}/usr/lib")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath=${SDK_NATIVE_SYSROOT}/lib")

# use alternative line below to support gdb debugger
#set(CMAKE_CXX_FLAGS "${ALLLANG_FLAGS} -ggdb -Og -ffunction-sections -Wall" CACHE STRING "" FORCE)
set(CMAKE_MAKE_PROGRAM "${SDK_NATIVE_SYSROOT}/usr/bin/make" CACHE FILEPATH "" FORCE)
# Qt
set(QT_HOST_PATH "${SDK_NATIVE_SYSROOT}")
set(QT_HOST_PATH_CMAKE_DIR "${SDK_NATIVE_SYSROOT}/usr/lib/cmake")
set(Python3_EXECUTABLE "${SDK_NATIVE_SYSROOT}/usr/bin/python3")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

set(CMAKE_NO_SYSTEM_FROM_IMPORTED true)
