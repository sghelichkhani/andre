//! \file
//! \brief Routines performing the actual IO
//!
//! This module contains the functions that perform the actual IO operations,
//! i.e. the ones that communicate with the operating system and the output
//! file. The module uses a sort of "attributes", i.e. static variables
//! globally accessible from everywhere within this module.


#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "vtkw_databox.h"
#include "vtkw_io.h"
#include "vtkw_abort.h"


//! Pointer to output file

//! We use a file pointer that is locally known within this module for
//! accessing the VTU file. File pointer must be initialised to NULL to allow
//! for testing if a file is currently open
static FILE *of = NULL;

//! Type of output file

//! In this variable we store the type of the output file, i.e. VTU or PVTU.
//! It is set by vtkwInitIO()
static vtkw_fileType ofType = NO_TYPE;


// ============
//  vtkwIOInit
// ============

//! Open output file and initisalise some XML tags

//! This function first attempts to open a file and then writes the file
//! header, i.e. the opening data independent XML tags of the VTK/Paraview
//! file. Currently it supports writing headers for .vtu and .pvtu files.
//!
//! \param fname        name of output file to open
//! \param ftype        type of vtk/paraview file
//! \param compressed   flag indicating whether data to be written was
//!                     compressed with zlib
void vtkwIOInit( char *fname, vtkw_fileType ftype, bool compressed ) {

  /* it is an error to open a new file before closing the old one */
  assert( of == NULL );

  /* open file */
  of = fopen( fname, "w" );
  if ( of == NULL ) {
    printf( "\n ERROR: Could not open '%s' for writing!\n\n", fname );
    vtkwAbort( NULL, VTKW_IO_ERROR );
  }

  /* write XML header */
  switch( ftype ) {

  case VTU:
    fprintf( of, "<?xml version=\"1.0\"?>\n" );
    fprintf( of, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" " );
    fprintf( of, "byte_order=\"LittleEndian\"" );
    if ( compressed ) {
      fprintf( of, " compressor=\"vtkZLibDataCompressor\"" );
    }
    fprintf( of, ">\n<UnstructuredGrid>\n" );
    break;

  case PVTU:
    fprintf( of, "<?xml version=\"1.0\"?>\n" );
    fprintf( of, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" " );
    fprintf( of, "byte_order=\"LittleEndian\">\n" );
    fprintf( of, "<PUnstructuredGrid GhostLevel=\"0\">\n" );
    break;

  case NO_TYPE:
  default:
    printf( "\n ERROR: Wrong type supplied to vtkwFileInit\n\n" );
    vtkwAbort( NULL, VTKW_IO_ERROR );
    break;
  }

  /* remember file type */
  ofType = ftype;

}


// ============
//  vtkwIOStop
// ============

//! Close some XML tags and output file

//! This function writes the file footer, i.e. it writes the closing tags of
//! those XML elements opened by vtkFileInit(). Currently it supports writing
//! the footers for .vtu and .pvtu files.
//! Additionally it closes the file and sets the file pointer to NULL.
void vtkwIOStop() {

  /* if we try to close something that's not open ... */
  assert( of != NULL );

  /* close XML tags */
  switch( ofType ) {

  case VTU:
    fprintf( of, "</UnstructuredGrid>\n" );
    fprintf( of, "</VTKFile>\n" );
    break;

  case PVTU:
    fprintf( of, "</PUnstructuredGrid>\n" );
    fprintf( of, "</VTKFile>\n" );
    break;

  case NO_TYPE:
  default:
    printf( "\n Severe Internal Error in %s at line %d!!!\n\n",
            __FILE__, __LINE__ );
    vtkwAbort( NULL, VTKW_IO_ERROR );
    break;
  }

  /* close file */
  if ( fclose(of) == EOF ) {
    printf( "\n ERROR: Problems closing VTU file!\n\n" );
    vtkwAbort( NULL, VTKW_IO_ERROR );
  }

  /* Nullify file pointer */
  of = NULL;

}


// ======================
//  vtkwIOWriteDataArray
// ======================

//! Output data array to file putting them inside a \<DataArray\> element

//! This function outputs the encoded and maybe also compressed data insides
//! a data box as a \<DataArray\> to the XMl file. Currently it supports only
//! writing .vtu files. The following table gives an overview on the machine
//! types specified by VTK that we use to represent the corresponding
//! conceptual types
//!
//! <center>
//! <table>
//! <tr><td><b>conceptual type (#vtkw_dataKind)</b></td><td><b>machine
//! type</b></td></tr>
//! <tr><td>#VTKW_POINTS           </td><td align="center"> Float32</td></tr>
//! <tr><td>#VTKW_CONNECTIVITY     </td><td align="center"> Int32  </td></tr>
//! <tr><td>#VTKW_OFFSETS          </td><td align="center"> Int32  </td></tr>
//! <tr><td>#VTKW_CELL_TYPES       </td><td align="center"> UInt8  </td></tr>
//! <tr><td>#VTKW_DIAMOND_ID       </td><td align="center"> UInt8  </td></tr>
//! <tr><td>#VTKW_PROCESS_ID       </td><td align="center"> UInt16 </td></tr>
//! <tr><td>#VTKW_POINT_DATA_SCALAR</td><td align="center"> Float32</td></tr>
//! <tr><td>#VTKW_POINT_DATA_VECTOR</td><td align="center"> Float32</td></tr>
//! <tr><td>#VTKW_CELL_DATA_SCALAR </td><td align="center"> Float32</td></tr>
//! <tr><td>#VTKW_CELL_DATA_VECTOR </td><td align="center"> Float32</td></tr>
//! </table>
//! </center>
//!
//! \param dBox   data box object containing the data
void vtkwIOWriteDataArray( const vtkw_dataBox *const dBox ) {

  Int32 nComponents = 0;
  Int32 nBytes = 0;
  char *vtkType = NULL;

  /* safety checks */
  assert( of != NULL );
  assert( ofType == VTU );
  assert( dBox->kind != VTKW_NO_KIND );

  /* determine data type and number of components */
  switch( dBox->kind ) {
  case VTKW_POINTS:
  case VTKW_POINT_DATA_VECTOR:
  case VTKW_CELL_DATA_VECTOR:
    nComponents = 3;
    vtkType = "Float32";
    break;
  case VTKW_CELL_TYPES:
    nComponents = 1;
    vtkType = "UInt8";
    break;
  case VTKW_CONNECTIVITY:
  case VTKW_OFFSETS:
    nComponents = 1;
    vtkType = "Int32";
    break;
  case VTKW_POINT_DATA_SCALAR:
  case VTKW_CELL_DATA_SCALAR:
    nComponents = 1;
    vtkType = "Float32";
    break;
  case VTKW_PROCESS_ID:
    nComponents = 1;
    vtkType = "UInt16";
    break;
  case VTKW_DIAMOND_ID:
    nComponents = 1;
    vtkType = "UInt8";
    break;
  default:
    printf( "\n Severe Internal Error in %s at line %d!!!\n\n",
            __FILE__, __LINE__ );
    vtkwAbort( NULL, VTKW_IO_ERROR );
    break;
  }

  /* open data array element */
  fprintf( of, "<DataArray type=\"%s\" Name=\"%s\" format=\"binary\"",
           vtkType, dBox->descr );
  fprintf( of, " NumberOfComponents=\"%1d\">\n", nComponents );

  /* output data */
  // for ( Int32 k = 0; k < dBox->eBufLen; k++ ) {
  //  fprintf( of, "%c", dBox->eBuf[k] );
  // }
  nBytes = fwrite( (void*)dBox->eBuf, 1, dBox->eBufLen, of );

  /* close data array element */
  fprintf( of, "</DataArray>\n" );
}


// ================
//  vtkwIOWriteTag
// ================

//! Low-level primitive for writing a tag to output file

//! This function is a low-level writing routine that simply inserts the
//! given string into the file followed by a line-break. Conceptually we
//! expect the string to be the opening or closing tag of an XML element.
//! However, it could also be a comment e.g.
//!
//! \param tag   string to be inserted
//!
//! \return nothing
void vtkwIOWriteTag( char *tag ) {

  assert( of != NULL );
  fprintf( of, "%s\n", tag );

}
