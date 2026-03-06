
compilation for CSD3

cmake ../. -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_BUILD_TYPE= -DCMAKE_SERIAL=False


## Diffraction patterns

Generate intensities at the diffracted beams using -opDP. 
Note that this is not a true diffraction pattern- as the intensity is only calculated at a set of predetermined g-vector points and not everywhere.
The intensity and size of spots are indicative only, these patterns can be used to check if we are in sensible diffraction conditions, rather than being intended for analytic work.

Four files are produced:

_before.diff.png        is the intensity at the g-vector positions at an idealised (defect free) foil back surface before rotating the foil.
_max.diff.png           is the maximum intensity anywhere in an idealised (defect free) . It is used in many beam mode to check where the energy is.
_ready.diff.png         is the intensity at the g-vector positions at an idealised (defect free) foil back surface after rotating the foil.
_after.diff.png         is the intensity at the g-vector positions actually used in the full calculation, using the atom positions used in the full calculation.


