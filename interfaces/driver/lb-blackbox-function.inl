/* lb-blackbox-function.inl
 * Copyright (C) 2017 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@lirmm.fr>
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


#include <lb-blackbox-function.h>
#include <lb-blackbox-abstract.h>
#include <lb-blackbox-functor.h>



/***********************************************************
 * API to launch a generic function over a linbox blackbox *
 **********************************************************/
extern BlackboxTable blackbox_hashtable;

// call a functor over a blackbox from the hashtable, result is given through 1st parameter
template<class Functor, class Result>
void BlackboxFunction::call(Result &res, const std::pair<const BlackboxKey, BlackboxAbstract*> &blackbox, const Functor &functor){
        ApplyBlackboxFunctor<Functor, Result> Ap(res, functor);
        (blackbox.second)->Accept(Ap);
}

// call a functor over a blackbox from the hashtable, no result
template<class Functor>
void BlackboxFunction::call(const std::pair<const BlackboxKey, BlackboxAbstract*> &v, const Functor &f){
        void *dumbresult;
        call(dumbresult,v,f);
}

// call a functor over a blackbox from its key, result is given through 1st parameter
template<class Functor, class Result>
void BlackboxFunction::call(Result &res, const BlackboxKey &key, const Functor &functor){
        BlackboxTable::const_iterator it = blackbox_hashtable.find(key);
        if (it != blackbox_hashtable.end())
                BlackboxFunction::call(res, *it, functor);
        else
                throw lb_runtime_error("LinBox ERROR: use of a non allocated blackbox \n");// throw an exception
}

// call a functor over a blackbox from its key, no result
template<class Functor>
void BlackboxFunction::call(const BlackboxKey &k, const Functor &f) {
        void *dumbresult;
        call(dumbresult,k,f);
}


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
