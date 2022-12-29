# Try to find OSQPLibrary

set(OSQP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/)

foreach(component
  libosqp
  libqdldl
)

  set(OSQP_LIBRARIES
    ${OSQP_LIBRARIES}
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/${component}.so)

endforeach()

if(OSQP_INCLUDE_DIR AND OSQP_LIBRARIES)

  set( OSQPLibrary_FOUND true )

endif()

if(OSQPLibrary_FOUND)
  message(STATUS "Found OSQP Library: ${OSQP_LIBRARIES}")
else()
  message(FATAL_ERROR "Could not find OSQP Library")
endif()
