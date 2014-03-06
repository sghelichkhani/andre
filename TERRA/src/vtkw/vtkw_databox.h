//! \file
//! \brief Most importantly definition of central #vtkw_dataBox structure.
//!
//! Apart from the usual prototes for functions from vtkw_databox.c this
//! file contains the defintion of the #vtkw_dataBox structure which is a
//! very central aspect of the VTKW library. This is accompanied by the
//! definition of the #dataBoxState enumeration type that is used for keeping
//! track of the internal state of the databox.

#ifdef VTKW_DATABOX_H
#else
#define VTKW_DATABOX_H

#include "vtkw_types.h"


//! Enumeration type for keeping track of the state of a vtkw_dataBox object

//! In order to ensure that the order of processing steps of data to be
//! exported is not violated we employ this enumeration data type to document
//! the current state of a vtkw_dataBox object. The state is updated each time
//! one of the following methods is called on the object.
typedef enum {

  //! no data inside; state after vtkwDataBoxCreate() or vtkwDataBoxFlush()
  DB_EMPTY,      
  DB_FILLED,     //!< data inside; state after vtkwDataBoxFill()
  DB_COMPRESSED, //!< compressed data inside; state after vtkwDataBoxCompress()
  DB_ENCODED,    //!< encoded data inside; state after vtkwDataBoxEncode()
  DB_WRITTEN     //!< data were written; state after vtkwDataBoxWrite()

} dataBoxState;

//! \brief Central data structure for VTKW library.
//!
//! The vtkw_databox type provides a container for managing the conversion,
//! compression, encoding and outputting of the actual data values that go
//! into the VTU file. Thus, it is probably the central data structure of
//! the VTKW implementation.
//!
//! The three steps
//! - conversion of field data / resp. generation of derived data
//! - compression of data
//! - encoding of data
//!
//! each require a target buffer for storing the result of the operation.
//! In order to reduce the memory footprint we only use two actual, i.e.
//! dynamically allocated, data buffers and associate the three virtual
//! target buffers #dBuf (data), #cBuf (compression) and #eBuf (encoding)
//! with them. The way this is done depends on the usage scenario, i.e.
//! whether we do compression or not:
//! 
//!    <table>
//!      <tr><td colspan="3">with compression:</td></tr>
//!        <tr><td>#dBuf</td><td>-></td><td>#bufB</td></tr>
//!        <tr><td>#cBuf</td><td>-></td><td>#bufA</td></tr>
//!        <tr><td>#eBuf</td><td>-></td><td>#bufB</td></tr>
//!    </table>
//!    and
//!    <table>
//!      <tr><td colspan="3">without compression:</td></tr>
//!        <tr><td>#dBuf</td><td>-></td><td>#bufA</td></tr>
//!        <tr><td>#cBuf</td><td>-></td><td>#bufA</td></tr>
//!        <tr><td>#eBuf</td><td>-></td><td>#bufB</td></tr>
//!    </table>
typedef struct {

  //! \brief Data buffer
  //!
  //! This is the buffer refered to as data buffer. It holds simulation and
  //! mesh data befored compression and encoding. It is a virtual buffer in
  //! the sense that it just points to #bufA or #bufB depending on the usage
  //! scenario.
  //!
  //! \note In case of field data the data are converted from double to single
  //!       precision when the databox is filled.
  unsigned char *dBuf;

  //! \brief Compression buffer
  //!
  //! This is the buffer refered to as compression buffer. It holds the data
  //! after they were compressed with zlib. It is a virtual buffer in the
  //! sense that it just points to #bufA or #bufB depending on the usage
  //! scenario.
  unsigned char *cBuf;

  //! \brief Encoding buffer
  //!
  //! This is the buffer refered to as encoding buffer. It holds the data
  //! after they were base64 encoded. It is a virtual buffer in the sense
  //! that it just points to #bufB.
  unsigned char *eBuf;

  //! \brief capacity of data buffer
  //!
  //! Capacity of data buffer, i.e. maximal number of entries that can be
  //! stored in it. Fixed when databox object is created.
  Int32 dBufCap;

  //! \brief capacity of compression buffer
  //!
  //! Capacity of compression buffer, i.e. maximal number of entries that can
  //! be stored in it. Fixed when databox object is created.
  Int32 cBufCap;

  //! \brief capacity of encoding buffer
  //!
  //! Capacity of encoding buffer, i.e. maximal number of entries that can
  //! be stored in it. Fixed when databox object is created.
  Int32 eBufCap;

  //! \brief current length of data buffer
  //!
  //! Length of data buffer, i.e. how many entries are currently stored in
  //! this buffer. May vary during lifetime of databox object.
  Int32 dBufLen;

  //! \brief current length of compression buffer
  //!
  //! Length of compression buffer, i.e. how many entries are currently stored
  //! in this buffer. May vary during lifetime of databox object.
  Int32 cBufLen;

  //! \brief current length of encoding buffer
  //!
  //! Length of encoding buffer, i.e. how many entries are currently stored
  //! in this buffer. May vary during lifetime of databox object.
  Int32 eBufLen;

  //! \brief Shorter of the two real auxilliary buffers
  //!
  //! We use two real buffers to represent the data, compression and encoding
  //! buffer. The two buffers need to have different sizes and this is the
  //! shorted one. Computation of the required maximal capacity and dynamic
  //! memory allocation is performed in the constructor metho
  //! vtkwDataBoxCreate(). Memory is deallocated in the destructor method
  //! vtkwDataBoxDestroy().
  unsigned char *bufA;

  //! \brief Longer of the two real auxilliary buffers
  //!
  //! We use two real buffers to represent the data, compression and encoding
  //! buffer. The two buffers need to have different sizes and this is the
  //! longer one. Computation of the required maximal capacity and dynamic
  //! memory allocation is performed in the constructor metho
  //! vtkwDataBoxCreate(). Memory is deallocated in the destructor method
  //! vtkwDataBoxDestroy().
  unsigned char *bufB;

  //! \brief Array with header info for \<DataArray\> element
  //!
  //! In the case that the data stored in a \<DataArray\> element are encoded
  //! and/or compressed we need to provide some additional information to VTK
  //! on the size of data before and after encoding/compression. We store this
  //! in this array. It also needs to be base64 encoded, so we do this in
  //! vtkwDataBoxEncode() and store it up front in the encoding buffer #eBuf.
  Int32 header[4];

  //! String describing data

  //! This is a string describing the data stored in the dataBox. It serves
  //! as the "Name" attribute of the \<DataArray\> element in the VTU file.
  //! In case of mandatory or derived data this description is hard-coded in
  //! the routine vtkwManagerExportData(). Apart from "diamondID" and
  //! "processID" the fixed values are mandated by the file format
  //! specification of the VTK. For the field data explicitely supplied
  //! by Terra we can freely chose a description in Terra. In total we have:
  //!
  //! <center>
  //! <table>
  //! <tr><td><b>type of data as #vtkw_dataKind</b></td><td><b>description
  //! string</b></td></tr>
  //! <tr><td>#VTKW_POINTS            </td><td> "points"       </td></tr>
  //! <tr><td>#VTKW_CONNECTIVITY      </td><td> "connectivity" </td></tr>
  //! <tr><td>#VTKW_OFFSETS           </td><td> "offsets"      </td></tr>
  //! <tr><td>#VTKW_CELL_TYPES        </td><td> "types"        </td></tr>
  //! <tr><td>#VTKW_DIAMOND_ID        </td><td> "diamondID"    </td></tr>
  //! <tr><td>#VTKW_PROCESS_ID        </td><td> "processID"    </td></tr>
  //! <tr><td>#VTKW_POINT_DATA_SCALAR </td><td> user's choice  </td></tr>
  //! <tr><td>#VTKW_POINT_DATA_VECTOR </td><td> user's choice  </td></tr>
  //! <tr><td>#VTKW_CELL_DATA_SCALAR  </td><td> user's choice  </td></tr>
  //! <tr><td>#VTKW_CELL_DATA_VECTOR  </td><td> user's choice  </td></tr>
  //! </table>
  //! </center>
  //! <br>
  char *descr;

  //! \brief Field to keep track of state of the databox object
  //!
  //! This field is used to keep track of the state of the databox object.
  //! The state is changed whenever a method, e.g. vtkwDataBoxFill() that
  //! changes or outputs the data stored in the object is executed. See
  //! description of #dataBoxState for further details.
  dataBoxState state;

  //! \brief Type of data stored inside databox
  //!
  //! This field describes the type of data stored inside the databox from
  //! the VTK perspective. See description of #vtkw_dataKind for further
  //! details.
  vtkw_dataKind kind;

} vtkw_dataBox;


/* Prototypes */
void vtkwDataBoxCreate( Int32 nNodes, Int32 nCells, bool compress );
void vtkwDataBoxDestroy();
void vtkwDataBoxFlush();
void vtkwDataBoxEncode();
void vtkwDataBoxCompress();
void vtkwDataBoxFill( const double *const data, vtkw_simSetup *simSet,
                      vtkw_dataKind dataType, char *descr );
void vtkwDataBoxWrite();

#endif
