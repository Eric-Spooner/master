
#-----------------------------------------------------------------------------
set(MODULE_NAME Logic)

#-----------------------------------------------------------------------------

#
# SlicerExecutionModel
#
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})



#-----------------------------------------------------------------------------
# Get ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
# Get VTK

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
  )

set(MODULE_SRCS
	benderIOUtils.cxx
	benderIOUtils.txx
	benderIOUtils.h
	benderWeightMap.cxx
	benderWeightMap.txx
	benderWeightMap.h
	benderWeightMapIO.h
	benderWeightMapIO.cxx
	benderWeightMapIO.txx
	benderWeightMapMath.h
  )

set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()
