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

/** @file matrix/blockedmatrix/dense-block-allocator.h
 * @ingroup blockedmatrix
 * A \c DenseBlockAllocator<_Field>\c allocated a uniform grid of
 * \cBlasMatrix<_Field>\c for a \c BlockedMatrix<_Field,_Allocator>\c.
 *
 */

#ifndef __LINBOX_matrix_blockedmatrix_dense_block_allocator_H
#define __LINBOX_matrix_blockedmatrix_dense_block_allocator_H

#include "linbox/matrix/blockedmatrix/block.h"
#include "linbox/matrix/dense-matrix.h"

#include <vector>

namespace LinBox
{
	/**
	 *
	 */
	template<class _Field>
	class DenseBlock : public Block<_Field>{
	public:
		typedef BlasMatrix<_Field> Repr;
	protected:
		Repr _repr;
	public:
		DenseBlock(const BlasMatrix<_Field>& repr) :
			_repr(repr){}

		virtual BlasMatrix<_Field>& operator*(){
			return _repr;
		}
	}; // end of class DenseBlock

	/**
	 *
	 */
	template<class _Field>
	class DenseBlockAllocator{
	public:
		typedef _Field              Field;
		typedef DenseBlock<Field> Block_t;

	protected:
		const Field*              _field;
		size_t                     _rows; //!< Number of rows of blocks
		size_t                     _cols; //!< Number of columns of blocks
		size_t               _block_rows; //!< Number of rows per block
		size_t               _block_cols; //!< Number of columns per block
		std::vector<Block_t> _store;

	public:
		//////////////////
		// CONSTRUCTORS //
		//////////////////

		/** Allocates \f$ M \times N\f$ blocks, with each block
		 * a new zero matrix \f$ m \times n\f$ in dimension.
		 * @param F
		 * @param m rows of blocks
		 * @param n cols of blocks
		 * @param block_m rows per block
		 * @param block_n cols per block
		 *
		 */
		DenseBlockAllocator(
			const Field& F,
			const size_t m,
			const size_t n,
			const size_t block_m,
			const size_t block_n) :
			// Init list begins here.
			_field(&F),
			_rows(m),
			_cols(n),
			_block_rows(block_m),
			_block_cols(block_n){}

		void allocate(){
			for(unsigned int i = 0; i < _rows * _cols; i++){
				_store.push_back(Block_t(BlasMatrix<Field>(*_field, _block_rows, _block_cols)));
			}

		}

		Block_t& lookupBlock(size_t i, size_t j){
			size_t offset = (_cols * i) + j;
			return _store[offset];
		}	

	}; // end of class DenseBlockAllocator
} // end of namespace LinBox

#endif // __LINBOX_matrix_blockedmatrix_dense_block_allocator_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
