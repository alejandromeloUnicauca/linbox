/*
 * Copyright (C) 2017 Matthew Wezowicz
 *
 * Written by :
 *               Matthew Wezowicz <mwezz@udel.edu>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/** @file matrix/blockedmatrix/block.h
 * @ingroup blockedmatrix
 * A \c Block<_Field>\c represents a single block of a 
 * \c BlockedMatrix<_Field,_Allocator>\C.
 *
 */

#ifndef __LINBOX_matrix_blockedmatrix_block_H
#define __LINBOX_matrix_blockedmatrix_block_H

namespace LinBox
{
	/**
	 *
	 */
	template<class _Field>
	class Block{
	public:
		virtual ~Block(){}
	}; // end of class Block
} // end of namespace LinBox

#endif // __LINBOX_matrix_blockedmatrix_block_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
