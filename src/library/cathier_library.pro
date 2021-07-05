TEMPLATE = lib
TARGET = guevara_cathier${BUILDMODEEXT}

INCBDIR =

#!include ../../config-cpp-library

HEADERS = \
# required for meshCleaner
   til/til_common.h \
   til/numeric_array.h \
   til/numeric_array_operators.h \
   til/numeric_array_tools.h \
   til/numeric_array.tpp \
   til/ext_declarations.h \
   til/labels.h \
   til/TExpr.h \
   til/TExprBasicFunctors.h \
   til/TExprConcatenation.h \
   til/TExprMacros.h \
   til/TExprPlaceHolders.h \
   til/TExprOperators.h \
   til/TExprFunctions.h \
   til/Matrix3.h \
   til/matrix3Tools.h \
   til/miscTools.h \
   til/misc_tools.tpp \
   til/templateTools.h \
   til/is_traits.h \
   til/cat2type.h \
   til/loop.h \
   til/basic_iterator.h \
   til/basic_range.h \
   til/functors.h \
   til/misc_scalar_functions.h \
   til/meta.h \
   til/multi_array.h \
   til/stditerator.h \
   til/std_wrap.h \
   til/traits.h \
   til/til_declarations.h \
   cathier/aims_wrap.h \
   cathier/misc_utils.tpp \
   cathier/misc_sort.h \
   cathier/misc_sort.tpp \
   cathier/mesh_decimation.h \
   cathier/globalTraits.h \
   cathier/miscUtils.h \
   cathier/mesh_conversion.h \
   cathier/Mesh.h \
   cathier/MeshTraits.h \
   cathier/meshUtils.h \
   cathier/mesh_conversion.tpp \
   cathier/mesh_utils.tpp \
   cathier/stditerator.h \
# required for constellation
   til/sparse_vector.h \
   til/sparse_vector_operators.h \
   til/sparse_vector_policies.h \
   til/sparse_vector_tools.h \
   til/value_proxy.h \
   til/value_proxy_policies.h \
   cathier/poly_solver.h \
   cathier/poly_solver_policies.h \
   cathier/poly_solver.tpp \
   cathier/geometrics.h \
   cathier/geometrics.tpp \
   cathier/kdtree.h \
   cathier/kdtree.tpp \
   cathier/triangle_mesh_geodesic_map.h \
   cathier/triangle_mesh_geodesic_map.tpp \
   cathier/fuzzy_logic.h \
   cathier/fuzzy_logic.tpp \
   cathier/declarations_external.h \
   cathier/SparseVector.h \
   cathier/graph_accessors.h \
   cathier/graph.h \
   cathier/fraction.h \
   cathier/fraction.tpp \
   cathier/fraction_policies.h \
   cathier/binary_tree.h \
   cathier/binary_tree.tpp \
   cathier/index_collection.h \
   cathier/cyclic_iterator.h

SOURCES = \
   cathier/aims_wrap.cpp

