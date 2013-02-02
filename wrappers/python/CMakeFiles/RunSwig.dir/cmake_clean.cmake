FILE(REMOVE_RECURSE
  "CMakeFiles/RunSwig"
  "../../python/src/swig_doxygen/OpenMMSwig.cxx"
  "../../python/simtk/openmm/openmm.py"
  "../../python/src/swig_doxygen/OpenMM.i"
  "../../python/src/swig_doxygen/swig_lib/python/pythonprepend.i"
  "../../python/src/swig_doxygen/swig_lib/python/pythonappend.i"
  "../../python/src/swig_doxygen/OpenMM_docstring.i"
  "../../python/src/swig_doxygen/OpenMM_headers.i"
  "../../python/src/swig_doxygen/swig_lib/python/exceptions.i"
  "../../python/src/swig_doxygen/swig_lib/python/extend.i"
  "../../python/src/swig_doxygen/swig_lib/python/features.i"
  "../../python/src/swig_doxygen/swig_lib/python/header.i"
  "../../python/src/swig_doxygen/swig_lib/python/pythoncode.i"
  "../../python/src/swig_doxygen/swig_lib/python/typemaps.i"
  "../../python/src/swig_doxygen/swig_lib/python/pythonprepend_all.i"
  "../../python/src/swig_doxygen/swigInputConfig.py"
  "../../python/src/swig_doxygen/swigInputBuilder.py"
  "src/swig_doxygen/doxygen/xml/index.xml"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/RunSwig.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
