cmake_minimum_required (VERSION 2.6)
project (KinFitter2)

# The version number.
set(KinFitter2_MAJOR_VERSION 1)
set(KinFitter2_MINOR_VERSION 0)
set(KinFitter2_PATCH_VERSION 0)
SET(KinFitter2_VERSION "${InSANE_MAJOR_VERSION}.${InSANE_MINOR_VERSION}.${InSANE_PATCH_VERSION}")

SET(lib_files
  KinFitter
  KinematicFitting
  )

#----------------------------------------------------------------------------
# Find Libraries
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist MathMore Graf Graf3d
   Gpad GX11 GeomPainter X3d Gviz3d RGL
   Eve EG Foam GenVector Geom GLEW Gui HistPainter MathCore Matrix MLP Net
   Physics Spectrum Thread Tree TreeViewer)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
#add_definitions(${ROOT_CXX_FLAGS})

find_package(PkgConfig REQUIRED)

list(APPEND CMAKE_PREFIX_PATH /group/clas12/packages/eigen/3.4.0/build)
pkg_search_module(Eigen3 REQUIRED eigen3)

#----------------------------------------------------------------------------
set(LINK_DIRECTORIES
   ${ROOT_LIBRARY_DIR}
   )
link_directories(${LINK_DIRECTORIES})

set(LINK_LIBRARIES 
   ${ROOT_LIBRARIES}
   )

 include_directories(
   ${Eigen3_INCLUDE_DIRS}
   ${ROOT_INCLUDE_DIR}
   ${CMAKE_BINARY_DIR}
   ${INCLUDE_DIRECTORIES}
   )

#@add_library(KinematicFitting cxx)

#----------------------------------------------------------------------------
set(aname              KinFitter2)
set(libname            "${aname}")
set(dictname           "${libname}Dict")
set(lib_LINKDEF        "${PROJECT_SOURCE_DIR}/include/LinkDef.h")
set(lib_DICTIONARY_SRC "${libname}Dict.cxx")
set(lib_PCM_FILE       "${PROJECT_BINARY_DIR}/lib${libname}Dict_rdict.pcm")

set(lib_SRCS)
set(lib_HEADERS)
foreach(infileName ${lib_files})
   SET(lib_SRCS    ${lib_SRCS}    ${PROJECT_SOURCE_DIR}/src/${infileName}.cxx)
   SET(lib_HEADERS ${lib_HEADERS} ${PROJECT_SOURCE_DIR}/include/${infileName}.h)
endforeach(infileName)
file(GLOB lib_HEADERS_HXX ${PROJECT_SOURCE_DIR}/include/*.hxx)
SET(lib_HEADERS ${lib_HEADERS} ${lib_HEADERS_HXX})

include_directories(${PROJECT_SOURCE_DIR}/include)

ROOT_GENERATE_DICTIONARY(${dictname} ${lib_HEADERS} LINKDEF ${lib_LINKDEF} OPTIONS -p)

SET(lib_HEADERS ${lib_HEADERS})
SET(lib_SRCS ${lib_Fortran_SRCs} ${lib_SRCS} ${lib_DICTIONARY_SRC})

SET(lib_VERSION "${${PROJECT_NAME}_VERSION}")
SET(lib_MAJOR_VERSION "${${PROJECT_NAME}_MAJOR_VERSION}")
SET(lib_LIBRARY_PROPERTIES 
    VERSION "${lib_VERSION}"
    SOVERSION "${lib_MAJOR_VERSION}"
    SUFFIX ".so")

ADD_CUSTOM_TARGET(${aname}_ROOTDICTS DEPENDS ${lib_SRCS} ${lib_HEADERS} ${lib_DICTIONARY_SRC} ${lib_DICTIONARY_HEADER})

add_library(${libname} SHARED ${lib_SRCS})
target_link_libraries(${libname} ${LINK_LIBRARIES} ${needs_libs})
set_target_properties(${libname} PROPERTIES ${lib_LIBRARY_PROPERTIES})
add_dependencies(${libname} ${needs_libs} ${aname}_ROOTDICTS)

install(
  TARGETS ${libname} 
  #   EXPORT ${PROJECT_NAME}Targets
  DESTINATION lib )

install(
  FILES ${lib_PCM_FILE} 
  DESTINATION lib )

install(
   FILES ${lib_HEADERS} 
   DESTINATION include/${PROJECT_NAME} )

