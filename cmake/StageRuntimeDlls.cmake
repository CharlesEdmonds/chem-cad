# Stages the full transitive DLL closure of one executable next to it.
#
# $<TARGET_RUNTIME_DLLS> only knows the DLLs CMake has imported targets for.
# Packaged RDKit pulls in Boost, cairo, freetype, zlib, expat and friends
# through plain PE imports that no CMake target describes, so an executable
# staged from that list still dies with STATUS_DLL_NOT_FOUND. Reading the
# import tables of the linked binary is the only complete answer.
#
# Several executables stage overlapping closures into the same directory in
# parallel, so each file is written under a unique temporary name and renamed
# into place; the loser of a race finds the file already there.
#
# Invoked as: cmake -D EXE=<binary> -D DEST=<dir> -D SEARCH=<dir> -P <this>

file(GET_RUNTIME_DEPENDENCIES
  EXECUTABLES "${EXE}"
  RESOLVED_DEPENDENCIES_VAR _resolved
  UNRESOLVED_DEPENDENCIES_VAR _unresolved
  CONFLICTING_DEPENDENCIES_PREFIX _conflict
  DIRECTORIES "${SEARCH}"
  PRE_EXCLUDE_REGEXES "api-ms-win-.*" "ext-ms-.*"
  POST_EXCLUDE_REGEXES ".*[/\\][Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\].*")

# A DLL already staged into DEST by an earlier target conflicts with the copy in
# SEARCH. They are the same file; always take the one from SEARCH so the size
# check below can short-circuit.
foreach(_name IN LISTS _conflict_FILENAMES)
  foreach(_candidate IN LISTS _conflict_${_name})
    string(FIND "${_candidate}" "${SEARCH}/" _fromSearch)
    if(_fromSearch EQUAL 0)
      list(APPEND _resolved "${_candidate}")
      break()
    endif()
  endforeach()
endforeach()

foreach(_dll IN LISTS _resolved)
  get_filename_component(_name "${_dll}" NAME)
  set(_target "${DEST}/${_name}")
  if(EXISTS "${_target}")
    file(SIZE "${_dll}" _sourceSize)
    file(SIZE "${_target}" _targetSize)
    if(_sourceSize EQUAL _targetSize)
      continue()
    endif()
  endif()

  string(RANDOM LENGTH 12 _suffix)
  set(_staged "${DEST}/${_name}.${_suffix}.tmp")
  file(COPY_FILE "${_dll}" "${_staged}" RESULT _copyError)
  if(_copyError)
    message(WARNING "chemcad: could not stage ${_name}: ${_copyError}")
    continue()
  endif()
  file(RENAME "${_staged}" "${_target}" RESULT _renameError)
  if(_renameError)
    # Another executable's staging step won the race; its copy is identical.
    file(REMOVE "${_staged}")
  endif()
endforeach()

if(_unresolved)
  list(REMOVE_DUPLICATES _unresolved)
  message(WARNING "chemcad: unresolved runtime dependencies for ${EXE}: ${_unresolved}")
endif()
