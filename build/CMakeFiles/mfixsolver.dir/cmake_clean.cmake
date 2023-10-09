file(REMOVE_RECURSE
  "mfixsolver"
  "mfixsolver.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/mfixsolver.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
