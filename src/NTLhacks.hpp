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
namespace NTL_NAMESPACE {
	void power(zz_p& x, const zz_p& a, const ZZ& e) {
		if (e == 0) x = 1;
		else if (a == 0 && e < 0)
			Error("zz_p: division by zero");
		else {
			long p = zz_p::modulus();
			long ee = e % (p-1);
			power(x, a, ee);
		}
	}
	void power(GF2& x, const GF2& a, const ZZ& e) {
		if (e == 0) x = 1;
		else if (a == 0 && e < 0)
			Error("GF2: division by zero");
		else x = a;
	}
}
