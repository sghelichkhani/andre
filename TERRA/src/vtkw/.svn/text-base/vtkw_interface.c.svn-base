//! \file
//! \brief Terra <-> VTKW interface routines
//!
//! This file contains routines for exporting data from Terra in XML format
//! following the VTK specification, so that the files can directly be read
//! by Paraview for post-processing. The routines in this module are the only
//! ones that are needed and should be called by Terra itself.
//!
//! In order to export data from Terra using VTKW the following order of
//! subroutine calls must be obeyed:
//!
//! <ol>
//! <li>Call VTKW_PREPARE(); this informs VTKW about essential simulation
//!     parameters like e.g. the meshing parameters and number of MPI
//!     processes</li>
//! <li>Call VTKW_FIELDS(); this triggers the actual writing of an output
//!     file containing the provided data</li>
//! </ol>
//!
//! This works in a single-shot fashion. In order to export another set of
//! data both functions must be called again and it is not possible to append
//! new data to an exisiting VTU file.
//!
//! \note
//! <ul>
//! <li>The Terra <-> VTKW interface does currently not follow the rules
//!     of C interoperability as were standardised with Fortran2003. The
//!     interface itself is still sort of experimental.</li>
//!
//! <li> The VTKW_FIELDS and VTKW_PREPARE macros in their current definitions
//!      work for compilation with GCC and ICC, other compilers may require
//!      different names for the C functions to be called from FORTRAN77</li>
//! </ul>


#include <stdbool.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include "vtkw_manager.h"


// These macros are used for compilation with GCC and ICC, other compilers may
// require different names for the C functions to be called from FORTRAN77
#define VTKW_FIELDS vtkw_fields_
#define VTKW_PREPARE vtkw_prepare_


// =============
//  VTKW_FIELDS
// =============

//! Output data from Terra to VTU file

//! Th is routine is used to export field data from Terra to a single VTU
//! file. The routine is variadic, i.e. the interface has no fixed number
//! of arguments. This allows us to export as many fields as we want e.g.
//! for quick debugging.
//!
//! \param bname          basename of the output files
//! \param numSclFields   number of nodal scalar fields to be written
//! \param numVecFields   number of nodal vector fields to be written
//! \param (...)          an arbitrary number of pointer pairs in the form
//!                       (double* field, char* marker); first all scalar
//!                       fields followed by all vector fields
void VTKW_FIELDS( char *bname, Int32 *numSclFields, Int32 *numVecFields, ... ){

  Int32 nsf = *numSclFields;
  Int32 nvf = *numVecFields;
  Int32 k;
  double **fieldList  = NULL;
  char   **markerList = NULL;

  /* get hold of fields and markers */
  va_list varParams;
  va_start( varParams, numVecFields );

  /* generate lists of fields and markers */
  fieldList  = calloc( nsf + nvf, sizeof(double*) );
  markerList = calloc( nsf + nvf, sizeof(char*) );
  for ( k = 0; k < nsf + nvf; k++ ) {
    fieldList [k] = va_arg( varParams, double* );
    markerList[k] = va_arg( varParams, char*   );
  }

  /* export fields */
  vtkwManagerExportData( fieldList, nsf, nvf, 0, 0, markerList, bname );

  /* cleanup */
  va_end( varParams );
  free( fieldList  );
  free( markerList );

}


// ==============
//  VTKW_PREPARE
// ==============

//! Prepare IO by passing information on simulation setup to VTKW library

//! This routine must be called before doing IO to provide the VTWK library
//! with required information on the simulation setup and the parallel
//! execution environment.
//!
//! \param numProcs      number of MPI processes
//! \param myNum         rank of MPI process executing this routine
//! \param nt            number of subdivisions along a subdomain edge
//!                      (spherical)
//! \param nr            number of subdivisions of the shell in radial
//!                      direction
//! \param nd            number of subdomains per process
//! \param nodes         array containing the cartesian coordinates of the
//!                      nodes
//!                      for the spherical discretisation of the subdomains
//!                      w.r.t. the unit sphere
//! \param radii         array containing the radii of the grid layers
//!                      (indexed inside in)
//! \param generatePVTU  flag signalling whether a .pvtu meta-file should be
//!                      generated or not
//! \param storeID       flag signalling whether for each cell of the mesh the
//!                      ID of the corresponding diamond and MPI process
//!                      should be inserted into the .vtu file
//! \param compress      flag signalling whether data should be compressed
//!                      with zlib
//!
//! \note The routine does not perform any actual work. Instead it passes
//!       all information onwards to vtkwManagerSetSimSetup().

void VTKW_PREPARE( Int32 *numProcs, Int32 *myNum, Int32 *nt, Int32 *nr,
                   Int32 *nd, double *nodes, double *radii,
                   Int32 *generatePVTU, Int32 *storeID, Int32 *compress ) {

  vtkwManagerSetSimSetup( *numProcs, *myNum, *nd, *nt, *nr, nodes, radii,
                          *generatePVTU != 0 ? true : false,
                          *storeID != 0 ? true : false,
                          *compress != 0 ? true : false );
}
