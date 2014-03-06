//! \file
//! \brief Attempt to cleanly abort in case of trouble
//!
//! While all actions of the VTKW library are performed purely sequentially
//! we must not forget that we are performing them within a parallel execution
//! environment. This module provides a function that attempts to cleanly
//! abort execution in the case of trouble.

#include "mpi.h"

// ===========
//  vtkwAbort
// ===========

//! Wrapper for MPI_Abort()

//! In a parallel execution environment we cannot simply exit() the program,
//! but should attempt to shutdown all MPI processes in a controlled manner.
//! Note that this is the only part of the VTKW library that actually depends
//! on the availability of an MPI library.
//!
//! \param mesg       error message (currently unused!)
//! \param errorCode  errorCode to return
void vtkwAbort( char *mesg, int errorCode ) {
  MPI_Abort( MPI_COMM_WORLD, errorCode );
}
