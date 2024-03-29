//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file superphoton.cpp
#include "superphoton.h"

/*!
 @param sp the superphoton to be copied
*/
superphoton & superphoton::operator = (const superphoton & sp) {
    if (&sp == this)
        return *this;
    
    this -> en = sp.en;
    this -> weight = sp.weight;
    this -> rnow = sp.rnow;
    this -> munow = sp.munow;
    this -> poldeg = sp.poldeg;
    this -> kmu = sp.kmu;
    this -> nscatter = sp.nscatter;
    this -> wp = sp.wp;
    
    return *this;
}
