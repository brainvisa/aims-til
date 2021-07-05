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
   cathier/aims_wrap.h \
   cathier/misc_utils.tpp \
   cathier/misc_sort.h \
   cathier/misc_sort.tpp \
   cathier/mesh_decimation.h \
   cathier/globalTraits.h \
   til/Matrix3.h \
   til/matrix3Tools.h \
# required for constellation
   til/sparse_vector.h \
   til/sparse_vector_operators.h \
   til/sparse_vector_policies.h \
   til/sparse_vector_tools.h \
   cathier/poly_solver.h \
   cathier/poly_solver_policies.h \
   cathier/poly_solver.tpp \
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
# unused ?
   til/Accumulator.h \
   til/basicFunctors.h \
   til/basic_iterator.h \
   til/basic_range.h \
   til/boost_wrap.h \
   til/cat2type.h \
   til/functors.h \
   til/ghost_array.h \
   til/if_then_else.h \
   til/is_traits.h \
   til/KeysInterpolation.h \
   til/loop.h \
   til/meta.h \
   til/misc_scalar_functions.h \
   til/miscTools.h \
   til/multi_array.h \
   til/NeighborhoodConfigurations.h \
   til/Neighborhood.h \
   til/neighborhoodTools.h \
   til/PointList.h \
   til/pointListTools.h \
   til/point_operators.h \
   til/point_tools.h \
   til/Ptr.h \
   til/range_tools.h \
   til/SmartObject.h \
   til/stditerator.h \
   til/std_wrap.h \
   til/stl_declarations.h \
   til/templateTools.h \
   til/til_declarations.h \
   til/til.h \
   til/til_math.h \
   til/traits.h \
   til/value_proxy.h \
   til/value_proxy_policies.h \
   til/misc_tools.tpp \
   cathier/accessors.h \
   cathier/binary_tree.h \
   cathier/cat2type.h \
   cathier/convert.h \
   cathier/cyclic_iterator.h \
   cathier/func_iterator.h \
   cathier/geometrics.h \
   cathier/index_collection.h \
   cathier/mesh_conversion.h \
   cathier/Mesh.h \
   cathier/MeshTraits.h \
   cathier/meshUtils.h \
   cathier/misc_sort.h \
   cathier/miscUtils.h \
   cathier/ordered_iterator.h \
   cathier/scalar_matrix.h \
   cathier/stditerator.h \
   cathier/std_wrap.h \
   cathier/binary_tree.tpp \
   cathier/mesh_conversion.tpp \
   cathier/mesh_utils.tpp

SOURCES = \
   cathier/aims_wrap.cpp

