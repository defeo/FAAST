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

#include "Types.hpp"
#include "Couveignes2000.hpp"
#include "FE-Liftup-Pushdown.hpp"
#include "FE-Trace-Frob.hpp"
#include "Field.hpp"
#include "FieldAlgorithms.hpp"
#include "FieldElement.hpp"
#include "GCD.hpp"
#include "FieldPolynomial.hpp"
#include "FieldPrecomputations.hpp"
#include "Minpols.hpp"
#include "utilities.hpp"
#include "NTLhacks.hpp"

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
