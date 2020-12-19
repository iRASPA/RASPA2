cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(simulate)

find_package(VTK COMPONENTS 
  vtkIOLegacy
  vtkCommonColor
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonMisc
  vtkFiltersGeneral
  vtkFiltersSources
  vtkFiltersHybrid
  vtkIOImage
  vtkInteractionStyle
  vtkRenderingCore
  vtkRenderingAnnotation
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkRenderingVolumeOpenGL2 QUIET)
if (NOT VTK_FOUND)
  message("Skipping simulate: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(simulate  main.cxx )
  target_link_libraries(simulate PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(simulate  main.cxx )
  target_link_libraries(simulate PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS simulate
    MODULES ${VTK_LIBRARIES}
    )
endif () 

