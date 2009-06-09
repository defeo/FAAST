/*
	This file is part of the FAAST library.

	Copyright (c) 2009 Luca De Feo and Éric Schost.

	The most recent version of FAAST is available at http://www.lix.polytechnique.fr/~defeo/FAAST

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; see file COPYING. If not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/
#include "faast.hpp"

#include "Types.c++"
#include "Couveignes2000.c++"
#include "FE-Liftup-Pushdown.c++"
#include "FE-Trace-Frob.c++"
#include "Field.c++"
#include "FieldAlgorithms.c++"
#include "FieldElement.c++"
#include "FieldPolynomial.c++"
#include "FieldPrecomputations.c++"
#include "Minpols.c++"
#include "utilities.c++"
#include "NTLhacks.c++"

using namespace FAAST;

/* Classes */
template class Field<zz_p_Algebra>;
template class FieldElement<zz_p_Algebra>;
template class FieldPolynomial<zz_p_Algebra>;

template class Field<ZZ_p_Algebra>;
template class FieldElement<ZZ_p_Algebra>;
template class FieldPolynomial<ZZ_p_Algebra>;

template class Field<GF2_Algebra>;
template class FieldElement<GF2_Algebra>;
template class FieldPolynomial<GF2_Algebra>;
