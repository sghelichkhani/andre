//! \file
//! \brief Only prototypes for functions in vtkw_manager.c
//!
//! This header file contains prototypes for functions in vtkw_manager.c,
//! but only for those that are "public" in the sense that other parts of
//! the library need to be able to call them.

#include "vtkw_types.h"

void vtkwManagerSetSimSetup( Int32 nProcs, Int32 myRank, Int32 nSubdoms,
                             Int32 nSphIntervals, Int32 nRadIntervals,
                             double *sphMesh, double *radii,
                             bool createPVTU, bool storeID, bool compress );

void vtkwManagerExportData( double **field,
                            Int32 numSclFieldsNodes, Int32 numVecFieldsNodes,
                            Int32 numSclFieldsCells, Int32 numVecFieldsCells,
                            char **marker, char *fname );

