set(sl_files shpopen.c safileio.c dbfopen.c)

set(sl_paths)
list(SORT sl_files)
list(REVERSE sl_files)
foreach(file IN ITEMS ${sl_files})
  #message(STATUS "Add ${file} to library program options")
  list(APPEND sl_paths ${sl_src}/${file})
endforeach()

set(lib my_shape_lib)
add_library(${lib} ${sl_paths})

