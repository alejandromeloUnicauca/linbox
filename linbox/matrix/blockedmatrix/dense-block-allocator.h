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
	// Forward declaration.
	template<class _Field>
	class DenseBlockAllocator;


	/**
	 *
	 */
	template<class _Field>
	class DenseBlock : public Block<_Field>{
	public:
		typedef _Field                  Field;
		typedef typename Field::Element Element;
		typedef BlasMatrix<Field>       Repr;

	protected:
		Repr* _repr;

	public:
		/**
		 *
		 */
		DenseBlock(
			BlockCoord coord,
			void* owner,
			Repr* repr) :
			// Init list begins here.
			Block<Field>::Block(BlockType::DENSE, coord, owner),
			_repr(repr){}
		/**
		 *
		 */
		DenseBlock(
			size_t i,
			size_t j,
			void* owner,
			Repr* repr) :
			// Init list begins here.
			Block<Field>::Block(BlockType::DENSE, i, j, owner),
			_repr(repr){}
		/**
		 *
		 */
		DenseBlock(const DenseBlock& rhs) :
			// Init list begins here.
			Block<Field>::Block(rhs),
			_repr(rhs._repr){}
		/**
		 *
		 */
		DenseBlock& operator=(const DenseBlock& rhs){
			_repr->copy(*(rhs._repr));
			return *this;
		}

		/**
		 *
		 */
		virtual ~DenseBlock(){}

		/**
		 *
		 */
		virtual Repr& operator*() const{
			return *_repr;
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
		typedef BlasMatrix<Field>    Repr;

	protected:
		const Field*               _field;
		size_t                      _rows; //!< Number of rows of blocks
		size_t                      _cols; //!< Number of columns of blocks
		size_t                _block_rows; //!< Number of rows per block
		size_t                _block_cols; //!< Number of columns per block
		std::vector<Block_t*> block_store;
		std::vector<Repr*>     repr_store;

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

		virtual ~DenseBlockAllocator(){
			deallocate();
		}

		void allocate(){
			for(size_t i = 0; i < _rows; i++){
				for(size_t j = 0; j < _cols; j++){
					Repr* matrix =
						new Repr(*_field,
							 _block_rows,
							 _block_cols);
					Block_t* block =
						new Block_t(i,
							    j,
							    this,
							    matrix);
					repr_store.push_back(matrix);
					block_store.push_back(block);
				}
			}
		}

		void deallocate(){
			for(size_t i = 0; i < _rows * _cols; i++){
				Repr* matrix = repr_store[i];
				Block_t* block = block_store[i];
				delete matrix;
				delete block;
			}
		}


		Block_t& lookupBlock(size_t i, size_t j){
			size_t offset = (_cols * i) + j;
			return *block_store[offset];
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
