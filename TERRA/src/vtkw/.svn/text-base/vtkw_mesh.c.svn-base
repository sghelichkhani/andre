//! \file
//! \brief "Low-level" functions that require specific Terra knowledge
//!
//! This module contains sort of "low-level" functions that require information
//! on the mesh discretisation of Terra and the current parallelisation
//! strategy, i.e. splitting into subdomains and number of MPI processes.
//! This information is required e.g. for computation of cartesian coordinates
//! of nodes or the computation of diamond IDs for cells.

#include <stdbool.h>
#include <assert.h>
#include "vtkw_types.h"
#include "vtkw_databox.h"


//! Macro mapping coordinate tuples to 1D address in Terra arrays.

//! Macro for determining the 1D index of the cc-th coordinate of the node
//! (i1,i2,id) in the spherical discretisation, or the 1D index of a single
//! grid node on a diamond on the cc-th layer, resp. the corresponding value
//! of a scalar field at that node.
#define IDXS(i1,i2,id,cc) ((nd*((cc)-1)+((id)-1))*(nt+1)+((i2)-1))*(nt+1)+(i1)

//! Macro mapping coordinate tuples to 1D address in Terra arrays.

//! Macro for determining the 1D index of the comp-th value of a vector field
//! at node (i1,i2,id,ir).
#define IDXV(i1,i2,id,ir,comp) (IDXS((i1),(i2),(id),(comp))+((ir)-1)*3*(nd)*(nt+1)*(nt+1))


// ======================
//  vtkwMeshSetCellTypes
// ======================

//! Write VTK cell type for wedge into buffer for each cell.

//! This function writes the cell types for the Terra mesh into a character
//! buffer. The cell type for Terra is constant, since each cell is treated
//! as a wedge. The corresponding cell type in VTK is number 13. We treat it
//! as a Uint8.
//!
//! \param buf      buffer into which to write the cell types
//! \param capacity capacity of write buffer
//! \param simSet   struct with meta information on simulation
//!
//! \return number of bytes written to buffer

Int32 vtkwMeshSetCellTypes( unsigned char *buf, Int32 capacity,
                            vtkw_simSetup *simSet ) {

  Int32 nt = simSet->nSphIntervals;
  Int32 nr = simSet->nRadIntervals;
  Int32 nd = simSet->nSubdoms;
  Int32 nCells = 2*nt*nt*nr*nd;

  assert( capacity >= nCells );

  const unsigned char cellType = 13;

  for ( Int32 k = 0; k < nCells; k++ ) {
    buf[k] = cellType;
  }

  return nCells;
}


// ====================
//  vtkwMeshSetOffsets
// ====================

//! Generate array with "offsets" for vtu file

//! This function generates the "offsets" data for the mesh description in the
//! vtu file. Since all cells have the same type, i.e. 13 which corresponds to
//! vtk_wedge, the offset is constantly 6. Thus, we just need to add up nCells
//! times that offset. The resulting values are written into the specified
//! output buffer. We use Uint32 as data type for this.
//!
//! \param buf      buffer into which to write the offsets
//! \param capacity capacity of write buffer
//! \param simSet   struct with meta information on simulation
//!
//! \return number of bytes written to buffer

Int32 vtkwMeshSetOffsets( unsigned char *buf, Int32 capacity,
                          vtkw_simSetup *simSet ) {

  Int32 nt = simSet->nSphIntervals;
  Int32 nr = simSet->nRadIntervals;
  Int32 nd = simSet->nSubdoms;
  Int32 nCells = 2*nt*nt*nr*nd;
  Int32 nBytes = 4*nCells;
  Int32 *outBuf = (Int32 *)buf;

  assert( capacity >= nBytes );

  const Int32 offset = 6;

  for ( Int32 k = 0; k < nCells; k++ ) {
    outBuf[k] = (k+1) * offset;
  }

  return nBytes;
}


// =========================
//  vtkwMeshSetConnectivity
// =========================

//! Writes connectivity information for cells into a buffer

//! This function computes the connectivity information for the cells of
//! Terra's and writes them into a character buffer. Connectity means which
//! nodes make up a cell. More precisely we need to specify for each cell
//! in a specified order the indices of the six nodes that make up this wedge.
//! For computing these indices we need to take into account that we output
//! the node coordinates just in their natural order, i.e. the way Terra
//! stores them in combination with Fortran's column-major array layout.
//! We use Int32 as data type for the indices. This is currently sufficient,
//! since in Terra we also use FORTRAN77's INTEGER data type.
//!
//! \param buf      buffer into which to write the connectivity
//! \param capacity capacity of write buffer
//! \param simSet   struct with meta information on simulation
//!
//! \return number of bytes written to buffer

Int32 vtkwMeshSetConnectivity( unsigned char *buf, Int32 capacity,
                               vtkw_simSetup *simSet ) {

  Int32 nt = simSet->nSphIntervals;
  Int32 nr = simSet->nRadIntervals;
  Int32 nd = simSet->nSubdoms;
  Int32 nCells = 2*nt*nt*nr*nd;
  Int32 nBytes = 24 * nCells;
  Int32 i1, i2, ir, id, k;
  Int32 *of = (Int32 *)buf;

  assert( capacity >= nBytes );

  k = 0;
  for ( ir = 1; ir <= nr; ir++ ) {
    for ( id = 1; id <= nd; id++ ) {
      for ( i2 = 1; i2 <= nt; i2++ ) {
        for ( i1 = 1; i1 <= nt; i1++ ) {

          /* upper prism/wedge associated with node */
          of[k++] = IDXS( i1  , i2  , id, ir   );
          of[k++] = IDXS( i1-1, i2+1, id, ir   );
          of[k++] = IDXS( i1-1, i2  , id, ir   );

          of[k++] = IDXS( i1  , i2  , id, ir+1 );
          of[k++] = IDXS( i1-1, i2+1, id, ir+1 );
          of[k++] = IDXS( i1-1, i2  , id, ir+1 );

          /* lower prism/wedge associated with node */
          of[k++] = IDXS( i1  , i2  , id, ir   );
          of[k++] = IDXS( i1  , i2+1, id, ir   );
          of[k++] = IDXS( i1-1, i2+1, id, ir   );

          of[k++] = IDXS( i1  , i2  , id, ir+1 );
          of[k++] = IDXS( i1  , i2+1, id, ir+1 );
          of[k++] = IDXS( i1-1, i2+1, id, ir+1 );
        }
      }
    }
  }

  return nBytes;
}


// ===================
//  vtkwMeshSetPoints
// ===================

//! Generate node coordinates and stores them in a buffer

//! This function computes the 3D cartesian coordinates of all nodes owned
//! by the current MPI process. This can easily be accomplished by using the
//! sphMesh and radii fields of the simSet input parameter. The coordinates
//! are converted to float and then stores in the designated buffer. The
//! nodes are processed in their natural ordering, i.e. in the order Terra
//! stores them in combination with Fortran's column-major array layout.
//!
//! \param buf      buffer into which to write the node coordinates
//! \param capacity capacity of write buffer
//! \param simSet   struct with meta information on simulation
//!
//! \return number of bytes written to buffer

Int32 vtkwMeshSetPoints( unsigned char *buf, Int32 capacity,
                         vtkw_simSetup *simSet ) {

  Int32 nt = simSet->nSphIntervals;
  Int32 nr = simSet->nRadIntervals;
  Int32 nd = simSet->nSubdoms;
  Int32 nNodes = (nt+1)*(nt+1)*(nr+1)*nd;
  Int32 nBytes = 12 * nNodes;
  Int32 i1, i2, ir, id, pos;

  double *nodes = simSet->sphMesh;
  double *radii = simSet->radii;
  float *dataBuf = (float *)buf;

  assert( capacity >= nBytes );

  pos = 0;
  for ( ir = 0; ir <= nr; ir++ ) {
    for ( id = 1; id <= nd; id++ ) {
      for ( i2 = 1; i2 <= nt+1; i2++ ) {
        for ( i1 = 0; i1 <= nt; i1++ ) {
          dataBuf[pos++] = (float)(radii[ir] * nodes[IDXS(i1,i2,id,1)]);
          dataBuf[pos++] = (float)(radii[ir] * nodes[IDXS(i1,i2,id,2)]);
          dataBuf[pos++] = (float)(radii[ir] * nodes[IDXS(i1,i2,id,3)]);
        }
      }
    }
  }
  
  return nBytes;
}


// ===========================
//  vtkwMeshSetScalarNodeData
// ===========================

//! Copy data from Terra array to buffer and convert to float on-the-fly

//! This function takes the values of a nodal scalar field generated in Terra,
//! converts them to float and writes the resulting values into the specified
//! buffer. Ordering of the data is left untouched.
//!
//! \param data      1D array with data of a nodal scalar field from Terra
//! \param buf       buffer into which to write the cell types
//! \param capacity  capacity of write buffer
//! \param simSet    struct with meta information on simulation
//!
//! \return number of bytes written to buffer

Int32 vtkwMeshSetScalarNodeData ( const double* const data, unsigned char *buf,
                                  Int32 capacity, vtkw_simSetup *simSet ) {

  float *dataBuf = (float*)buf;
  Int32 nt = simSet->nSphIntervals;
  Int32 nr = simSet->nRadIntervals;
  Int32 nd = simSet->nSubdoms;
  Int32 nNodes = (nt+1)*(nt+1)*(nr+1)*nd;
  Int32 nBytes = 4 * nNodes;
  // Int32 i1, i2, id, ir
  Int32 pos;

  assert( capacity >= nBytes );

  // pos = 0;
  // for ( ir = 1; ir <= nr+1; ir++ ) {
  //   for ( id = 1; id <= nd; id++ ) {
  //     for ( i2 = 1; i2 <= nt+1; i2++ ) {
  //       for ( i1 = 0; i1 <= nt; i1++ ) {
  //         dataBuf[pos++] = (float)data[IDXS(i1,i2,id,ir)];
  //       }
  //     }
  //   }
  // }

  for ( pos = 0; pos < nNodes; pos++ ) {
    dataBuf[pos] = (float)data[pos];
  }

  return nBytes;
}


// ===========================
//  vtkwMeshSetVectorNodeData
// ===========================

//! Copy data from Terra array to buffer and convert to float on-the-fly

//! This function takes the values of a nodal scalar field generated in Terra,
//! converts them to float and writes the resulting values into the specified
//! buffer. Ordering of the data is left untouched.
//!
//! \param data     1D array with data of a nodal vector field from Terra
//! \param buf      buffer into which to write the cell types
//! \param capacity capacity of write buffer
//! \param simSet   struct with meta information on simulation
//!
//! \return   number of bytes written to buffer

Int32 vtkwMeshSetVectorNodeData ( const double* const data, unsigned char *buf,
                                  Int32 capacity, vtkw_simSetup *simSet ) {

  float *dataBuf = (float*)buf;
  Int32 nt = simSet->nSphIntervals;
  Int32 nr = simSet->nRadIntervals;
  Int32 nd = simSet->nSubdoms;
  Int32 nNodes = (nt+1)*(nt+1)*(nr+1)*nd;
  Int32 nBytes = 12 * nNodes;
  Int32 i1, i2, id, ir, pos, cc;

  assert( capacity >= nBytes );

  pos = 0;
  for ( ir = 1; ir <= nr+1; ir++ ) {
    for ( id = 1; id <= nd; id++ ) {
      for ( i2 = 1; i2 <= nt+1; i2++ ) {
        for ( i1 = 0; i1 <= nt; i1++ ) {
          for ( cc = 1; cc <= 3; cc++ ) {
            dataBuf[pos++] = (float)data[IDXV(i1,i2,id,ir,cc)];
          }
        }
      }
    }
  }

  return nBytes;
}


// ======================
//  vtkwMeshSetProcessID
// ======================

//! Determine and store MPI process ranks for cells

//! This function can be used to fill the databox with information on the
//! partitioning of the mesh into subdomains. We generate a scalar field
//! associated to the cells and simply set for each porocess its process id.
//!
//! We make the assumption that the process id can be represented by a 16-bit
//! unsigned integral number. This is checked by an assertion if NDEBUG macro
//! is not set.
//!
//! \param buf       buffer into which to write the cell process IDs
//! \param capacity  capacity of write buffer
//! \param simSet    struct with meta information on simulation
//!
//! \return number of bytes written to buffer

Int32 vtkwMeshSetProcessID( unsigned char *buf, Int32 capacity,
                            vtkw_simSetup *simSet ) {

  Int32 nt = simSet->nSphIntervals;
  Int32 nr = simSet->nRadIntervals;
  Int32 nd = simSet->nSubdoms;
  Int32 nCells = 2*nt*nt*nr*nd;
  Int32 nBytes = 2 * nCells;

  Uint16 *of = (Uint16 *)buf;

  assert( UINT16_MAX >= simSet->myRank );
  assert( capacity >= nBytes );

  for ( Int32 k = 0; k < nBytes; k++ ) {
    of[k] = (Uint16)simSet->myRank;
  }

  return nBytes;
}


// ======================
//  vtkwMeshSetDiamondID
// ======================

//! Determine and store diamond IDs for cells

//! This function can be used to fill the databox with a cell-based scalar
//! field that contains for each cell the number of the diamond it belongs
//! to. The diamond numbers are elements of [1,10], so we can use Uint8 for
//! their representation.
//!
//! \param buf        buffer into which to write the diamond IDs
//! \param capacity   capacity of write buffer
//! \param simSet     struct with meta information on simulation
//!
//! \return number of bytes written to buffer

Int32 vtkwMeshSetDiamondID( unsigned char *buf, Int32 capacity,
                            vtkw_simSetup *simSet ) {

  Int32 nt = simSet->nSphIntervals;
  Int32 nr = simSet->nRadIntervals;
  Int32 nd = simSet->nSubdoms;
  Int32 nCells = 2*nt*nt*nr*nd;
  Int32 nBytes = nCells;
  Int32 ii, id, ir, pos;
  unsigned char offset = 0;

  assert( capacity >= nBytes );

  if ( nd == 5 && simSet->myRank >= (simSet->nProcs/2) ) {
    offset = 5;
  }

  pos = 0;
  for ( ir = 1; ir <= nr; ir++ ) {
    for ( id = 1 + offset; id <= nd + offset; id++ ) {
      for ( ii = 1; ii <= 2*nt*nt; ii++ ) {
        buf[pos++] = (unsigned char)id;
      }
    }
  }

  return nBytes;
}
