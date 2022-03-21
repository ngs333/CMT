#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "edlib::edlib" for configuration ""
set_property(TARGET edlib::edlib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(edlib::edlib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib64/libedlib.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS edlib::edlib )
list(APPEND _IMPORT_CHECK_FILES_FOR_edlib::edlib "${_IMPORT_PREFIX}/lib64/libedlib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
