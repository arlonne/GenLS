version: v1.1

GenLS (Generate Lexicographical Structures for materials design)
is a tool for building models in Lexicographical Order!

Pro to Random generation:
1) search all permutation without repeatation and missing error
2) symmetries can be added

Pre to Random generation:
1) Current, it required many disk to store the full-permutation in large system.
However, this could be improved in the future.

In the future:
1) mpi/omp parallel version.
2) split and joint the file Allpermutation.out. 
The permutaion can be start from any old order. This is useful in large systems.
3) symmetries to decrease the number of permutations in defects mode.
4) working on interface with abinit codes (VASP,ATK and so on) and start calculations automatically.
5) more structure opteration should be added, such as rotation, move on a vector, and so on.
