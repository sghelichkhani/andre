//! \file
//! \brief Implementation of methods for central #vtkw_dataBox object
//!
//! The #vtkw_dataBox object provides a container for managing the conversion,
//! compression, encoding and exporting of data to the output VTU file. In this
//! module we define as "class" attribute a pointer to the databox object and
//! implement its methods. Conceptually the databox methods are one level
//! below the vtkwManager "class" in our implementation hierarchy.

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>

#include "vtkw_databox.h"
#include "vtkw_encode.h"
#include "vtkw_types.h"
#include "vtkw_abort.h"
#include "vtkw_mesh.h"
#include "vtkw_io.h"


//! Pointer to databox object

//! This variable is a pointer to the central databox object. The pointer
//! is set by the vtkwDataBoxCreate() constructor method and nullified by
//! the vtkwDataBoxDestroy() destructor method. Via this pointer the object
//! and its attributes are directly accessible by all methods defined in this
//! module.
//! \note For performance and implementation reasons we break encapsulation
//!       a little by passing the pointer to the object itself or the buffers
//!       it contains to selected methods on the lower levels of our conceptual
//!       hierarchy.
static vtkw_dataBox* dBox = NULL;

// internal prototype
char* vtkwDataSetDescr( const char *const srcStr );


//! \name Public Member functions:
//! These functions need to be visible to other parts of the VTKW library.
//! Thus, their prototypes are specified in the header file vtkw_databox.h
//@{

// ===================
//  vtkwDataBoxCreate
// ===================

//! Constructor method for #vtkw_dataBox object

//! This function creates a #vtkw_dataBox object and initialises it to a
//! pristine state, i.e. the struct and its auxilliary arrays are dynamically
//! allocated and the remaining parameters set to an initial state. After
//! construction the internal state of the object is #DB_EMPTY. The generated
//! object is accessible via the #dBox pointer from everywhere in this "class".
//!
//! \param nNodes    total number of nodes for which MPI process is responsible
//! \param nCells    total number of cells for which MPI process is responsible
//! \param compress  flag indicating whether compression of data is to be
//!                  performed before encoding and outputting them 
void vtkwDataBoxCreate( Int32 nNodes, Int32 nCells, bool compress ) {

  Int32 bufSizeA, bufSizeB;
  uLong auxVal, maxDatLen, maxAlloc;

  assert( dBox == NULL );

  // allocate single databox object
  dBox = (vtkw_dataBox*) malloc( sizeof(vtkw_dataBox) );
  if ( dBox == NULL ) {
    printf( "\n ERROR: Failed to alloc databox object (nNodes=%d)!\n\n",
            nNodes );
    vtkwAbort( NULL, VTKW_DATABOX_ERROR );
  }
 
  // initialise first set of components
  dBox->bufA      = NULL;
  dBox->bufB      = NULL;
  dBox->dBuf      = NULL;
  dBox->cBuf      = NULL;
  dBox->eBuf      = NULL;
  dBox->dBufLen   = 0;
  dBox->cBufLen   = 0;
  dBox->eBufLen   = 0;
  dBox->descr     = NULL;
  dBox->state     = DB_EMPTY;
  dBox->kind      = VTKW_NO_KIND;
  dBox->header[0] = 0;
  dBox->header[1] = 0;
  dBox->header[2] = 0;
  dBox->header[3] = 0;


  /**********************************************************************
   *  Determine maximal required length for the two auxilliary buffers  *
   **********************************************************************/

  // Item 1:
  // We are pessimistic and assume that compression is either not performed
  // or does not reduce the data size. Thus, the sizes are directly related
  // to the maximal amount of data per node resp. cell.
  //
  // Item 2:
  // We also need to take into account that base64 encoding enlarges the
  // size required by a factor of 4/3.
  //
  // Item 3:
  // Using compress from the zlib library, the destination buffer, upon entry,
  // must be at least 0.1% larger than length of the source data plus 12 bytes.
  // We can check this requirement dynamically by calling compressBound()
  maxDatLen = 12 * nNodes > 24 * nCells ? 12 * nNodes : 24 * nCells;
  auxVal = compressBound( maxDatLen );
  maxAlloc = ( auxVal > maxDatLen ) ? (Int32)auxVal : (Int32)maxDatLen;
  bufSizeA = (Int32)maxAlloc;
  assert( maxAlloc == (uLong)bufSizeA );
  bufSizeB = ( (bufSizeA + 2) / 3) * 4 + 24;
  assert( bufSizeA > 0 && bufSizeB > 0 );

  // allocate the shorter of the two auxilliary buffers
  dBox->bufA = (unsigned char*) malloc( bufSizeA );
  if ( dBox->bufA == NULL ) {
    printf( "\n ERROR: Failed to alloc (%f MB) memory for databox\n\n",
            (float)bufSizeA / 1024.0 / 1024.0 );
    vtkwAbort( NULL, VTKW_DATABOX_ERROR );
  }

  // allocate the longer of the two auxilliary buffers
  dBox->bufB = (unsigned char*) malloc( bufSizeB );
  if ( dBox->bufB == NULL ) {
    printf( "\n ERROR: Failed to alloc (%f MB) memory for databox\n\n",
            (float)bufSizeB / 1024.0 / 1024.0 );
    vtkwAbort( NULL, VTKW_DATABOX_ERROR );
  }

  // depending on usage (compress or not) we have different associations
  // of the conceptual buffers to the memory buffers
  if ( compress ) {
    dBox->dBuf = dBox->bufB;
    dBox->cBuf = dBox->bufA;
    dBox->eBuf = dBox->bufB;

    dBox->dBufCap = bufSizeB;
    dBox->cBufCap = bufSizeA;
    dBox->eBufCap = bufSizeB;
  }
  else {
    dBox->dBuf = dBox->bufA;
    dBox->cBuf = dBox->bufA;
    dBox->eBuf = dBox->bufB;

    dBox->dBufCap = bufSizeA;
    dBox->cBufCap = bufSizeA;
    dBox->eBufCap = bufSizeB;
  }
}


// ====================
//  vtkwDataBoxDestroy
// ====================

//! Destructor method for #vtkw_dataBox object

//! This function destroys a #vtkw_dataBox object and frees its auxilliary
//! arrays and string components. It re-sets the #dBox pointer to NULL.
void vtkwDataBoxDestroy() {

  assert( dBox != NULL );
  free( dBox->bufA );
  free( dBox->bufB );
  free( dBox->descr );
  free( dBox );
  dBox = NULL;

}


// ==================
//  vtkwDataBoxFlush
// ==================

//! Flush contents of databox

//! This function can be used to discard the contents of the databox object.
//! It should only be called after the contents were written. Note that we
//! do not de-allocate memory. Thus, the capacity remains unchanged.
//! Calling this method changes the state of the object to #DB_EMPTY.
void vtkwDataBoxFlush() {

  assert( dBox != NULL && dBox->state == DB_WRITTEN );
  dBox->dBufLen = 0;
  dBox->eBufLen = 0;
  dBox->cBufLen = 0;
  dBox->kind   = VTKW_NO_KIND;
  dBox->state  = DB_EMPTY;

}


// ===================
//  vtkwDataBoxEncode
// ===================

//! Triggers encoding of data stored in the databox

//! This function triggers base64 encoding of the data stored in the data
//! box object. The data to be encoded is taken from the #vtkw_dataBox::cBuf
//! array and the result written to the #vtkw_dataBox::eBuf array. Calling
//! this method changes the state of the object to #DB_ENCODED.
void vtkwDataBoxEncode() {

  unsigned char *hdr = NULL;

  assert( dBox->state == DB_FILLED || dBox->state == DB_COMPRESSED );

  // Case 1: with compression
  if ( dBox->state == DB_COMPRESSED ) {

    // First encode the header
    hdr = (unsigned char*)(dBox->header);
    vtkwEncodeBase64Buffer2Buffer( hdr, 16, dBox->eBuf );

    // Now encode data */
    vtkwEncodeBase64Buffer2Buffer( dBox->cBuf, dBox->cBufLen, dBox->eBuf+24 );

    // Set new data length
    dBox->eBufLen = (dBox->cBufLen+2)/3*4 + 24;
  }

  // Case 2: no compression
  else {

    // Data reside in dBuf
    assert( dBox->dBuf == dBox->cBuf );
    dBox->cBufLen = dBox->dBufLen;
    dBox->header[2] = dBox->dBufLen;
    dBox->header[3] = dBox->dBufLen;

    // First encode the header
    hdr = (unsigned char*)(dBox->header+2);
    vtkwEncodeBase64Buffer2Buffer( hdr, 4, dBox->eBuf );

    // Now encode data */
    vtkwEncodeBase64Buffer2Buffer( dBox->cBuf, dBox->cBufLen, dBox->eBuf+8 );

    // Set new data length
    dBox->eBufLen = (dBox->cBufLen+2)/3*4 + 8;
  }

  dBox->state = DB_ENCODED;

}


// =====================
//  vtkwDataBoxCompress
// =====================

//! Perform compression of data stored in the databox

//! This function performs compression of the data stored in the data box
//! object by means of the zlib library. The data to be encoded is taken from
//! the #vtkw_dataBox::dBuf array and the result written to the
//! #vtkw_dataBox::cBuf array. Calling this method changes the state of the
//! object to #DB_COMPRESSED.
void vtkwDataBoxCompress() {

  uLong dstLen = (uLong)dBox->cBufCap;
  uLong srcLen = (uLong)dBox->dBufLen;
  int retCode;

  /* Perform compression */
  assert( dBox->state == DB_FILLED );
  assert( dBox->cBuf != dBox->dBuf );
  assert( compressBound( srcLen ) <= dstLen );
  retCode = compress( dBox->cBuf, &dstLen, dBox->dBuf, srcLen );
  assert( retCode == Z_OK );

  // Update information
  //
  // We compressed all data in one go. Thus, from a vtk perspective we have
  // 1 block and blocksize is srcLen. Since there is only one block, general
  // blocksize and blocksize of last block are identical. We set the header
  // correspondingly.
  dBox->cBufLen   = dstLen;
  dBox->header[0] = 1;
  dBox->header[1] = srcLen;
  dBox->header[2] = srcLen;
  dBox->header[3] = dstLen;
  dBox->state     = DB_COMPRESSED;
}


// ==================
//  vtkwDataBoxWrite
// ==================

//! Triggers export of data to output file

//! This function can be used to output the encoded and maybe compressed data
//! stored in the databox object. Actual work is delegated to the function
//! vtkwIOWriteDataArray(). Calling this method changes the state of the
//! object to #DB_WRITTEN.
void vtkwDataBoxWrite() {

  assert( dBox != NULL && dBox->state == DB_ENCODED );
  vtkwIOWriteDataArray( dBox );
  dBox->state = DB_WRITTEN;

}


// ==================
//  vtkwDataBoxFill
// ==================

//! Insert field data into databox

//! This function takes a given scalar or vector field and inserts its values
//! into the databox object. In doing so the data a converted from double to
//! float and re-interpreted as char in turn. The function should only be used
//! for double precision data.
//!
//! We assume that the data buffer contains data from Terra in Fortran layout.
//! Thus, we use a mapping macro to find the corresponding 1D index.
//!
//! Calling this method changes the state of the object to #DB_FILLED.
//!
//! \param data      pointer to data array containg values of the field
//! \param simSet    struct with meta data of the simulation run
//! \param dataType  type of data to be converted
//! \param descr     textual identifier of data
void vtkwDataBoxFill( const double *const data, vtkw_simSetup *simSet,
                      vtkw_dataKind dataType, char *descr ) {

  Int32 nt = simSet->nSphIntervals;
  Int32 nr = simSet->nRadIntervals;
  Int32 nd = simSet->nSubdoms;
  Int32 nNodes = (nt+1)*(nt+1)*(nr+1)*nd;

  /* Perform some safety checks */
  assert( descr != NULL );
  assert( sizeof(float) == 4 );
  assert( dBox->state == DB_EMPTY );
  assert( dBox->dBufCap >= 12*nNodes );

  /* We need to distinguish between nodal and cell data and
     scalar and vector fields */
  switch( dataType ) {

  case VTKW_POINT_DATA_SCALAR:
    dBox->dBufLen = vtkwMeshSetScalarNodeData( data, dBox->dBuf, dBox->dBufCap,
                                               simSet );
    break;

  case VTKW_POINT_DATA_VECTOR:
    dBox->dBufLen = vtkwMeshSetVectorNodeData( data, dBox->dBuf, dBox->dBufCap,
                                               simSet );
    break;

  case VTKW_POINTS:
    dBox->dBufLen = vtkwMeshSetPoints( dBox->dBuf, dBox->dBufCap, simSet );
    break;

  case VTKW_CELL_TYPES:
    dBox->dBufLen = vtkwMeshSetCellTypes( dBox->dBuf, dBox->dBufCap, simSet );
    break;

  case VTKW_OFFSETS:
    dBox->dBufLen = vtkwMeshSetOffsets( dBox->dBuf, dBox->dBufCap, simSet );
    break;

  case VTKW_CONNECTIVITY:
    dBox->dBufLen = vtkwMeshSetConnectivity( dBox->dBuf, dBox->dBufCap,
                                             simSet );
    break;

  case VTKW_PROCESS_ID:
    dBox->dBufLen = vtkwMeshSetProcessID( dBox->dBuf, dBox->dBufCap, simSet );
    break;

  case VTKW_DIAMOND_ID:
    dBox->dBufLen = vtkwMeshSetDiamondID( dBox->dBuf, dBox->dBufCap, simSet );
    break;

  case VTKW_CELL_DATA_SCALAR:
  case VTKW_CELL_DATA_VECTOR:
    printf( "\n ERROR: Transformation of cell data not implemented, yet!\n\n");
    vtkwAbort( NULL, VTKW_DATABOX_ERROR );
    break;

  default:
    printf( "\n ERROR: vtkwDataBoxFill() called with dataType = %d !!\n\n",
            dataType );
    vtkwAbort( NULL, VTKW_DATABOX_ERROR );
  }

  /* Set remaining information */
  dBox->descr = vtkwDataSetDescr( descr );
  dBox->kind  = dataType;
  dBox->state = DB_FILLED;

}

//@}


//! \name Private Member functions:
//! These functions are only used by the vtkw_datbox internally. They are
//! hidden from the outside by not putting their prototypes into the header
//! file vtkw_databox.h
//@{

// ==================
//  vtkwDataSetDescr
// ==================

//! Auxilliary function for string business

//! Small auxilliary function for handling the business of inserting the
//! description of the data into the databox object. Some of the strings
//! describing the conceptual type of the data inside the databox (see
//! #vtkw_dataBox::descr for details) are hard-coded string constants.
//! We dynamically re-alloc the dBox->descr string and perform an sprintf()
//! operation in order to avoid problems with this.
//!
//! \param srcStr string to be written into description
char* vtkwDataSetDescr( const char *const srcStr ) {

  char *newStr = NULL;
  newStr = (char*) realloc( dBox->descr, strlen(srcStr) + 1 );
  sprintf( newStr, "%s", srcStr );

  return newStr;
}
