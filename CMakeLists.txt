
#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(doiPET)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Setup of ROOT ANALYSIS : optional. 
option(WITH_ANALYSIS_USE "Build example with analysis objects" OFF)
if(WITH_ANALYSIS_USE)
add_definitions(-DUSEROOT)
else()
endif()

option(DICOM_USE_DCMTK "Build with DCMTK support" ON)
if(DICOM_USE_DCMTK)
  # Set DCMTK root for Homebrew
  set(DCMTK_ROOT "/opt/homebrew/opt/dcmtk")
  set(DCMTK_INCLUDE_DIRS "${DCMTK_ROOT}/include")
  
  # Add library directory to linker search path
  link_directories("${DCMTK_ROOT}/lib")
  
  # Manually specify DCMTK libraries with full paths
  set(DCMTK_LIBRARIES 
      "${DCMTK_ROOT}/lib/libdcmrt.dylib"
      "${DCMTK_ROOT}/lib/libdcmimgle.dylib"
      "${DCMTK_ROOT}/lib/libdcmimage.dylib"
      "${DCMTK_ROOT}/lib/libdcmdata.dylib"
      "${DCMTK_ROOT}/lib/libdcmiod.dylib"
      "${DCMTK_ROOT}/lib/libdcmfg.dylib"
      "${DCMTK_ROOT}/lib/liboflog.dylib"
      "${DCMTK_ROOT}/lib/libofstd.dylib"
      "${DCMTK_ROOT}/lib/liboficonv.dylib"
  )
  
  add_definitions(-DDICOM_USE_DCMTK)
  add_definitions(-DG4_DCMTK)
  include_directories(${DCMTK_INCLUDE_DIRS})
  message(STATUS "DCMTK found: ${DCMTK_INCLUDE_DIRS}")
  message(STATUS "DCMTK libraries: ${DCMTK_LIBRARIES}")
endif()
#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)
if(WITH_ANALYSIS_USE)
	EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
endif(WITH_ANALYSIS_USE) 
#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/shared/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
                  ${PROJECT_SOURCE_DIR}/shared/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
                  ${PROJECT_SOURCE_DIR}/shared/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(doiPET doiPET.cc ${sources} ${headers})

if(WITH_ANALYSIS_USE)
	EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LD_FLAGS})
endif(WITH_ANALYSIS_USE)

target_link_libraries(doiPET ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
if(DICOM_USE_DCMTK)
  target_link_libraries(doiPET ${DCMTK_LIBRARIES})
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(DOIPET_SCRIPTS
  init_vis.mac
  vis.mac
  run_imageQualityPhantom_wholeBody.mac
  run_imageQualityPhantom_smallAnimal.mac
  run_NECR.mac
  run_spatialResolution.mac
  run_sensitivity.mac 
  run_normalization.mac 
  run.mac
  inputParameter.txt
  )

foreach(_script ${DOIPET_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

install(TARGETS doiPET DESTINATION bin )
