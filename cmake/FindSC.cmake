#.rst
# FindSC
# ---------
#
# Find the sc library.
#
# This module sets the following variables:
#
# ::
#
#   SC_FOUND - Set to true if the library has been found.
#   SC_LIBRARIES - The libraries needed to use sc.
#   SC_INCLUDE_DIRS - The include directories needed to use sc.
#
# To provide the module with a hint to where the sc library has been
# installed, you can set the SC_ROOT environment variable to the required
# location.


# Copyright Maison de la Simulation (CEA/CNRS/INRIA/Univ. Paris Sud/UVSQ) - USR3441 CNRS
# Copyright EM2C (Ecole Centrale Paris) - UPR288 CNRS
#
# CanoP is a versatile software package designed for solving
# computational fluid dynamics problems using a cell-based
# adaptive mesh refinement approach.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

include (FindPackageHandleStandardArgs)

find_path (SC_INCLUDE_DIR
        NAMES sc.h
        HINTS ENV SC_ROOT C_INCLUDE_PATH)
mark_as_advanced (SC_INCLUDE_DIR)

find_library (SC_LIBRARY
        NAMES sc
        HINTS ENV SC_ROOT LD_LIBRARY_PATH)
mark_as_advanced (SC_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set SC_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SC DEFAULT_MSG SC_LIBRARY SC_INCLUDE_DIR)

if (SC_FOUND)
    find_package (ZLIB REQUIRED)
    find_package (Lua51 REQUIRED)

    set (SC_LIBRARIES ${SC_LIBRARY} ${ZLIB_LIBRARIES} ${LUA_LIBRARIES})
    set (SC_INCLUDE_DIRS ${SC_INCLUDE_DIR} ${ZLIB_INCLUDE_DIRS})

    #
    #   Only add MPI libs if sc was compiled with MPI support
    #
    file (STRINGS
            "${SC_INCLUDE_DIR}/sc_config.h"
            SC_WITH_MPI
            REGEX "^#define SC*_MPI")

    if (SC_WITH_MPI)
        find_package (MPI REQUIRED)
        set (SC_LIBRARIES ${SC_LIBRARIES} ${MPI_LIBRARIES})
        set (SC_INCLUDE_DIRS ${SC_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
    endif (SC_WITH_MPI)

    #
    #   Only add LAPACK libs if sc was compiled with LAPACK support
    #
    file (STRINGS
            "${SC_INCLUDE_DIR}/sc_config.h"
            SC_WITH_LAPACK
            REGEX "^#define SC*_LAPACK")

    if (SC_WITH_LAPACK)
        find_package (LAPACK REQUIRED)
        set (SC_LIBRARIES ${SC_LIBRARIES} ${LAPACK_LIBRARIES})
    endif (SC_WITH_LAPACK)
endif (SC_FOUND)
