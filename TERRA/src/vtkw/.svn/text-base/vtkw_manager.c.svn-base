//! \file
//! \brief High-level routines for data manipulation and export
//!
//! In our hierarchical layout the vtkwManager "class" is on the second level
//! just below the vtkwInterface "class" that constitutes the actual interface
//! between Terra and the VTKW library. As can nicely be seen e.g. from the
//! call graph of the vtkwManagerExportData() method the vtkwManager "class"
//! acts as a high-level driver and steers e.g.
//!
//! - generation of the vtkwDataBox object
//! - concersion, compression and encoding of data
//! - computation of derived data such as e.g. cell connectivities
//! - generation of filenames
//! - writing to the output file
//!
//! It does this by relying on methods of the underlying classes on deeper
//! levels of the hierarchy.

#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include "vtkw_types.h"
#include "vtkw_databox.h"
#include "vtkw_io.h"


//! Attribute storing information on simulation run

//! This attribute of type vtkw_simSetup is used to store information on the
//! current simulation run of Terra. It is initialised to unsensible default
//! values for error checking and filled with real information by the
//! vtkwManagerSetSimSetup() method.
static vtkw_simSetup simSet = { -1, -1, -1, -1, -1, NULL, NULL,
                                false, false, false };


//! \name Private Member functions:
//! These functions are only used by the vtkw_manager internally. They are
//! hidden from the outside by not putting their prototypes into the header
//! file vtkw_manager.h
//@{


// ==============================
//  vtkwManagerGenerateFileNames
// ==============================

//! Generate filenames for VTU and PVTU files

//! The function generates filenames for the VTU, and if necessary, also for
//! the PVTU file. The filenames are generated from the given basename argument
//! in combination with the information stored in the #simSet attribute. The
//! filenames are generated following the rule
//!
//! <center><b><code>\<baseName\>-\<myRank\>.vtu</code></b></center>
//! <center><b><code>\<baseName\>.pvtu</code></b></center>
//!
//! where <dfn>\<myRank\></dfn> is a field of length four containing the MPI
//! rank of the current process with leading zeros, if neccessary.
//!
//! \param baseName        basename of output files
//! \param vtuFileName     name of .vtu file for this MPI process
//! \param pvtuFileName    name of meta-file, only generated for process #0
//!                        and, if feature was enabled
//!
//! \note Memory to hold the filenames is dynamically allocated. It is the
//!       caller's responsibility to clear it, once no longer needed.
//!       If this process is not #0 or the feature was not enabled, the
//!       pvtuFileName parameter is not touched!
void vtkwManagerGenerateFileNames( char *baseName, char **vtuFileName,
                                   char **pvtuFileName ) {

  // Generate name of vtu-file for this MPI process
  *vtuFileName = (char*)malloc( strlen(baseName)+10 );
  sprintf( *vtuFileName, "%s-%4.4d.vtu", baseName, simSet.myRank );

  // Generate name of meta-file for process 0
  if ( simSet.createPVTU && simSet.myRank == 0 ) {
    *pvtuFileName = (char*)malloc( strlen(baseName)+6 );
    sprintf( *pvtuFileName, "%s.pvtu", baseName );
  }
}


// =======================
//  vtkwManagerHandleData
// =======================

//! Driver routine for processing steps and output of data in databox

//! This function accumulates all calls to methods of the databox object
//! required to output data in one place. It acts as a primitive driver
//! routine performing in that order the following steps
//!
//! <ol>
//! <li>Compression of data by using vtkwDataBoxCompress()</li>
//! <li>Encoding of data by using vtkwDataBoxEncode()</li>
//! <li>Writing of data to output file by using vtkwDataBoxWrite()</li>
//! <li>Flushing of databox object by using vtkwDataBoxFlush()</li>
//! </ol>
//!
//! \param data      pointer to data array containg values of the field
//! \param dataType  type of data to be converted
//! \param descr     textual identifier of data
void vtkwManagerHandleData( const double *const data, vtkw_dataKind dataType,
                            char *descr ) {

  vtkwDataBoxFill( data, &simSet, dataType, descr );
  if ( simSet.compress ) {
    vtkwDataBoxCompress();
  }
  vtkwDataBoxEncode();
  vtkwDataBoxWrite();
  vtkwDataBoxFlush();

}
//@}



//! \name Public Member functions:
//! These functions need to be visible to other parts of the VTKW library.
//! Thus, their prototypes are specified in the header file vtkw_manager.h
//@{


// ========================
//  vtkwManagerSetSimSetup
// ========================

//! Fill #simSet attribute with basic information on Terra run

//! The function inserts the central mesh and subdomain parameters and values
//! of the current simulation run into the #simSet attribute.
//!
//! \param nProcs          total number of MPI processes
//! \param myRank          MPI rank of process executing function
//! \param nSubdoms        number of subdomains associated to MPI process
//! \param nSphIntervals   number of subdivisions of a subdomain edge
//!                        (spherical)
//! \param nRadIntervals   number of subdivisions of the shell in radial
//!                        direction
//! \param sphMesh         1D array containing cartesian coordinates of the
//!                        nodes owned by the process w.r.t. to Terra's
//!                        discretisation of the unit sphere
//! \param radii           1D array containing the positions of the radial
//!                        layers in Terra's mesh
//! \param createPVTU      flag indicating whether a .pvtu parallel meta-file
//!                        should also be written for Paraview
//! \param storeID         flag toggling output of a unique diamond and MPI
//!                        process ID for each mesh cell
//! \param compress        flag indicating whether compression of data is
//!                        desired or not
void vtkwManagerSetSimSetup( Int32 nProcs, Int32 myRank, Int32 nSubdoms,
                             Int32 nSphIntervals, Int32 nRadIntervals,
                             double *sphMesh, double *radii,
                             bool createPVTU, bool storeID, bool compress ) {

  simSet.nProcs        = nProcs;
  simSet.myRank        = myRank;
  simSet.nSubdoms      = nSubdoms;
  simSet.nSphIntervals = nSphIntervals;
  simSet.nRadIntervals = nRadIntervals;
  simSet.sphMesh       = sphMesh;
  simSet.radii         = radii;
  simSet.createPVTU    = createPVTU;
  simSet.storeID       = storeID;
  simSet.compress      = compress;
}


// =======================
//  vtkwManagerExportData
// =======================

//! High-level driver for generation of VTU file and export of data

//! This function acts as a high-level driver for the generation of the VTU
//! file and the export of the explicitely (fields) or implicitely (nodal
//! coordinates, ...) given data. It is responsible for triggering of derived
//! data and the processing, i.e. conversion, compression and encoding, of
//! provided data fields and the export of those to the outpout file.
//!
//! \param field              2D array containing pointers to the fields to be
//!                           exported; fields must be in the following order:
//!                           - nodal scalar fields
//!                           - nodal vector fields
//!                           - cellular scalar fields
//!                           - cellular vector fields
//! \param numSclFieldsNodes  number of nodal scalar fields in field array
//! \param numVecFieldsNodes  number of nodal vector fields in field array
//! \param numSclFieldsCells  number of cellular scalar fields in field array
//! \param numVecFieldsCells  number of cellular vector fields in field array
//! \param marker             2D array with pointers to strings specifying the
//!                           names of the fields in the same order as in field
//!                           array
//! \param bname              basename of output files
void vtkwManagerExportData( double **field,
                            Int32 numSclFieldsNodes, Int32 numVecFieldsNodes,
                            Int32 numSclFieldsCells, Int32 numVecFieldsCells,
                            char **marker, char *bname ) {

  char tag[256];
  Int32 nNodesPerSubdom, nCellsPerSubdom, nNodes, nCells;
  Int32 curField = 0, k;
  char *vtuFileName  = NULL;
  char *pvtuFileName = NULL;


  /******************
   *  Preparations  *
   ******************/

  // Compute derived values
  nNodesPerSubdom = (simSet.nSphIntervals+1) * (simSet.nSphIntervals+1)
    * (simSet.nRadIntervals+1);
  nNodes = nNodesPerSubdom * simSet.nSubdoms;
  nCellsPerSubdom = 2 * simSet.nSphIntervals * simSet.nSphIntervals
    * simSet.nRadIntervals;
  nCells = nCellsPerSubdom * simSet.nSubdoms;

  // Make sure that vtkwManagerSetSimSetup() was called
  assert( simSet.nProcs > 0 );

  // Generate a data box for handling the individual fields
  vtkwDataBoxCreate( nNodes, nCells, simSet.compress );

  // Generate filenames
  vtkwManagerGenerateFileNames( bname, &vtuFileName, &pvtuFileName );

  // Open .vtu file and write XML header stuff
  vtkwIOInit( vtuFileName, VTU, simSet.compress );  

  // Write opening tag for piece element
  sprintf( tag, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">",
           nNodes, nCells );
  vtkwIOWriteTag( tag );


  /*********************
   *  Write Mesh Info  *
   *********************/

  // Points
  vtkwIOWriteTag( "<Points>" );
  vtkwManagerHandleData( NULL, VTKW_POINTS, "points" );
  vtkwIOWriteTag( "</Points>" );

  // Cells
  vtkwIOWriteTag( "<Cells>" );
  vtkwManagerHandleData( NULL, VTKW_CONNECTIVITY, "connectivity" );
  vtkwManagerHandleData( NULL, VTKW_OFFSETS, "offsets" );
  vtkwManagerHandleData( NULL, VTKW_CELL_TYPES, "types" );
  vtkwIOWriteTag( "</Cells>" );


  /********************
   *  Write PointData *
   ********************/
  if ( numSclFieldsNodes > 0 || numVecFieldsNodes ) {

    // open PointData element */
    vtkwIOWriteTag( "<PointData>" );

    // loop over scalar fields
    for ( k = 0; k < numSclFieldsNodes; k++ ) {
      vtkwManagerHandleData( field[curField], VTKW_POINT_DATA_SCALAR,
                             marker[curField] );
      curField++;
    }

    // loop over vector fields
    for ( k = 0; k < numVecFieldsNodes; k++ ) {
      vtkwManagerHandleData( field[curField], VTKW_POINT_DATA_VECTOR,
                             marker[curField] );
      curField++;
    }

    // close PointData element
    vtkwIOWriteTag( "</PointData>" );
  }

  
  /********************
   *  Write CellData  *
   ********************/
  if ( numSclFieldsCells > 0 || numVecFieldsCells || simSet.storeID ) {

    // open point data element
    vtkwIOWriteTag( "<CellData>" );

    // if user requires, we write also data on process and diamond ids
    if ( simSet.storeID ) {
      vtkwManagerHandleData( NULL, VTKW_DIAMOND_ID, "diamondID" );
      vtkwManagerHandleData( NULL, VTKW_PROCESS_ID, "processID" );
    }

    // close point data element
    vtkwIOWriteTag( "</CellData>" );
  }


  /*********************
   *  Finish vtu-File  *
   *********************/

  // close piece element
  vtkwIOWriteTag( "</Piece>" );
 
  // finalise XML tree structure and close file
  vtkwIOStop();


  /***************
   *  Meta File  *
   ***************/

  // We make process 0 responsible for generating the .pvtu meta file
  if ( simSet.createPVTU && simSet.myRank == 0 ) {
    vtkwIOInit( pvtuFileName, PVTU, simSet.compress );

    // Info on point data
    vtkwIOWriteTag( "<PPointData>" );

    curField = 0;
    for ( k = 0; k < numSclFieldsNodes; k++ ) {
      sprintf( tag, "<PDataArray type=\"Float32\" Name=\"%s\"/>",
               marker[curField++] );
      vtkwIOWriteTag( tag );
    }

    for ( k = 0; k < numVecFieldsNodes; k++ ) {
      sprintf( tag, "<PDataArray type=\"Float32\" Name=\"%s\"",
               marker[curField++] );
      sprintf( tag + strlen(tag), " NumberOfComponents=\"3\"/>" );
      vtkwIOWriteTag( tag );
    }

    vtkwIOWriteTag( "</PPointData>" );

    // Info on cell data
    vtkwIOWriteTag( "<PCellData>" );
    if ( simSet.storeID ) {
      sprintf( tag, "<PDataArray type=\"UInt8\" Name=\"diamondID\"" );
      sprintf( tag + strlen(tag), " NumberOfComponents=\"1\"/>" );
      vtkwIOWriteTag( tag );
      sprintf( tag, "<PDataArray type=\"UInt16\" Name=\"processID\"" );
      sprintf( tag + strlen(tag), " NumberOfComponents=\"1\"/>" );
      vtkwIOWriteTag( tag );
    }
    vtkwIOWriteTag( "</PCellData>" );

    // Info on points of the mesh
    vtkwIOWriteTag( "<PPoints>" );
    sprintf( tag, "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>");
    vtkwIOWriteTag( tag );
    vtkwIOWriteTag( "</PPoints>" );

    // Where to find the pieces
    for ( k = 0; k < simSet.nProcs; k++ ) {
      sprintf( tag, "<Piece Source=\"%s-%4.4d.vtu\"/>", bname, k );
      vtkwIOWriteTag( tag );
    }

    vtkwIOStop( );  
  }


  /*************
   *  Cleanup  *
   *************/
  vtkwDataBoxDestroy();
  free(  vtuFileName );
  free( pvtuFileName );

}

//@}
