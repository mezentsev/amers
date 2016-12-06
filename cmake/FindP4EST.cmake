#.rst
# FindP4EST
# ---------
#
# Find the p4est library.
#
# This module sets the following variables:
#
# ::
#
#   P4EST_FOUND - Set to true if the library has been found.
#   P4EST_LIBRARIES - The libraries needed to use p4est.
#   P4EST_INCLUDE_DIRS - The include directories needed to use p4est.
#
# To provide the module with a hint to where the p4est library has been
# installed, you can set the P4EST_ROOT environment variable to the required
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

find_path (P4EST_INCLUDE_DIR
        NAMES p4est.h
        HINTS ENV P4EST_ROOT C_INCLUDE_DIR)
mark_as_advanced (P4EST_INCLUDE_DIR)

find_library (P4EST_LIBRARY
        NAMES p4est
        HINTS ENV P4EST_ROOT LD_LIBRARY_PATH)
mark_as_advanced (P4EST_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set P4EST_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args (P4EST DEFAULT_MSG
        P4EST_LIBRARY
        P4EST_INCLUDE_DIR)

if (P4EST_FOUND)
    find_package(SC REQUIRED)

    set (P4EST_LIBRARIES ${P4EST_LIBRARY}
            ${SC_LIBRARIES})
    set (P4EST_INCLUDE_DIRS ${P4EST_INCLUDE_DIR}
            ${SC_INCLUDE_DIRS})
endif (P4EST_FOUND)