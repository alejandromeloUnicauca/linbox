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

#include <utility>

#ifndef __LINBOX_matrix_blockedmatrix_block_H
#define __LINBOX_matrix_blockedmatrix_block_H

namespace LinBox
{
	/**
	 *
	 */
	enum BlockType{
		DENSE
	};

	/**
	 *
	 */
	class BlockCoord{
	protected:
		size_t row;
		size_t col;

		BlockCoord();
	public:
		/**
		 *
		 */
		BlockCoord(const size_t i, const size_t j) :
			// Init list begins here.
			row(i),
			col(j){}

		/**
		 *
		 */
		BlockCoord(const BlockCoord& rhs) :
			// Init list begins here.
			row(rhs.row),
			col(rhs.col){}

		/**
		 *
		 */
		BlockCoord& operator=(const BlockCoord& rhs){
			row = rhs.row;
			col = rhs.col;
			return *this;
		}

		/**
		 *
		 */
		virtual ~BlockCoord(){}

		/**
		 *
		 */
		const size_t getRow(){
			return row;
		}

		/**
		 *
		 */
		const size_t getCol(){
			return col;
		}

		/**
		 *
		 */
		std::pair<size_t, size_t> getIJ(){
			return std::pair<size_t, size_t>(row, col);
		}

		/**
		 *
		 */
		bool operator==(const BlockCoord& rhs){
			return (row == rhs.row) && (col == rhs.col);
		}

		/**
		 *
		 */
		bool operator!=(const BlockCoord& rhs){
			return !(*this == rhs);
		}
	}; // end of class BlockCoord

	/**
	 *
	 */
	template<class _Field>
	class Block{
	public:
		typedef _Field Field;
		typedef typename Field::Element Element;

	protected:
		BlockType blockType;
		BlockCoord blockCoord;
		void* _owner;

		Block();
	public:
		/**
		 *
		 */
		Block(BlockType type, BlockCoord coord, void* owner) :
			// Init list begins here.
			blockType(type),
			blockCoord(coord),
			_owner(owner){}

		/**
		 *
		 */
		Block(BlockType type, size_t i, size_t j, void* owner) :
			// Init list begins here.
			blockType(type),
			blockCoord(i,j),
			_owner(owner){}

		/**
		 *
		 */
		Block(const Block& rhs) :
			// Init list begins here.
			blockType(rhs.blockType),
			blockCoord(rhs.blockCoord),
			_owner(rhs.owner){}

		/**
		 *
		 */
		Block& operator=(const Block& rhs);

		/**
		 *
		 */
		virtual ~Block(){}

		/**
		 *
		 */
		const BlockType getBlockType(){
			return blockType;
		}

		/**
		 *
		 */
		const BlockCoord& getBlockCoord(){
			return blockCoord;
		}

		/**
		 *
		 */
		bool operator==(const Block& rhs){
			bool same = blockType == rhs.blockType;
			same &= blockCoord == rhs.blockCoord;
			same &= _owner == rhs._owner;
			return same;
		}

		/**
		 *
		 */
		bool operator!=(const Block& rhs){
			return !(*this == rhs);
		}
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
