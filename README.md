Copyright (C) 2009-2013 Luca De Feo and Eric Schost.

FAAST is an open-source C++ library providing strongly object-oriented
data structures and high-performance algorithms for elements and
polynomials in Artin-Schreier towers over finite fields.

The latest stable version of FAAST is available at
https://github.com/defeo/FAAST/archive/latest.zip.

You can grab the latest sources from the GitHub repo at
https://github.com/defeo/FAAST.

The API documentation is available at https://defeo.github.io/faast,
or can be downloaded from
https://github.com/defeo/FAAST/archive/gh-pages.zip.


Copying Conditions:
===================

FAAST is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.

FAAST is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with FAAST; see the file COPYING. If not, write to the Free
Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.


Installing Instructions:
========================

Requirements
------------

FAAST is built on top of Victor Shoup's NTL, you can download it at
http://www.shoup.net/NTL/. Only NTL 5.x is supported, work to port the
library to NTL 6 is ongoing.

There's a few libraries that affect how NTL works. GMP is a
multiprecision arithmetic library that can be used by NTL to represent
large integers. Even though it is not expected to affect the overall
performances of FAAST, compiling NTL with GMP is still the recommended
method. GMP can be downloaded at http://gmplib.org/.

gf2x is a library for multiplying polynomials over the binary
field. If you work on characteristic 2 fields and you are interested
in performance, then we recommend that you install gf2x and compile
NTL with support for it. gf2x can be downloaded at
http://gforge.inria.fr/projects/gf2x/.

Instructions on how to compile NTL with support for GMP and/or gf2x
can be found at http://www.shoup.net/ntl/doc/tour.html.

Compling
--------

To compile FAAST, simply type

	tar xzf faast-xxx.tgz
	cd faast-xxx/
	./configure
	make
	make check

During the configure script FAAST checks for the presence of NTL,
GMP and gf2x; checks for GMP and gf2x are optional and can be disabled
through `--disable-ntl-static`. If NTL has been built as a shared library,
passing `--disable-ntl-static` might speed up the compilation.

Notice that you might have to specify the paths to NTL, GMP 
or gf2x if they are installed in some exotic directory; see `configure --help`.

Installing
----------

To install FAAST, from the compilation directory type

	su
	make install

The library will be installed to a default location, you can
change this by giving a different `--prefix` to the configure script.

Generating the API docs
-----------------------

The latest API docs are available online at
https://defeo.github.io/faast, or can be downloaded from
https://github.com/defeo/FAAST/archive/gh-pages.zip.

If you wish, you can generate from your source pacakge. You must have
installed the doxygen documentation system for this. You can download
doxygen at http://www.doxygen.org/. You will also need LaTeX and
Xy-pic to generate some formulae, graphs and two logos.

Once doxygen is installed, you can generate the documentation by typing

	make doc

This creates a directory doc/html containing the documentation. You
can edit the file `doxy.conf` to enable more options and output formats,
refer to the doxygen manual.

If you are interested in understanding the internals of FAAST, you can
generate a more detailed documentation by typing

	make doc-dev

This creates a directory doc-dev/html containing the documentation.
