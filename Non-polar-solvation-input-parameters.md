---
layout: page
---
        

#### Input Keywords for Non-polar Solvation Energy

##### apolar

    apolar            = yes


This will specify whether to do non-polar calculation. Value `yes` will allow the calculation wheras value `no` will not allow the calculation.


#### <a name="sasa_model"></a> SASA model

##### gamma

    gamma            =  0.02267

This specifies the surface tension proportionality term (kJ mol <sup>-1</sup> Å <sup>-2</sup> ) of the solvent. With `gamma = 0`, SASA energy is not calculated.

***

##### sasaconst

    sasaconst        = 3.84928

Offset or constant `c` in (kJ mol<sup>-1</sup>) from _E_ =  _γ_ * _SASA_ + _c_

***

##### sasrad

    sasrad           = 1.4

Solvent probe radius in Å to calculate solvent accessible surface area.

***

#### <a name="sav_model"></a> SAV model

***

##### press

    press            = 0.234304

This specifies the solvent pressure proportionality term (kJ mol<sup>-1</sup> Å<sup>-3</sup>). With `press = 0`, SAV energy is not calculated.

***

##### savrad

    savrad           = 1.29

Solvent probe radius in Å to calculate solvent accessible volume.

***

##### savconst

    savconst         = 0

Offset or constant `c` in kJ/mol from _E_ = ( <var>p</var> * _SAV_ ) + <var>c</var>

***

#### <a name="wca_like"> </a> Continuum/Integral based model (WCA-like)

##### WCA

    WCA              = yes

Switch for WCA method: "yes" or "no"

***

##### wcarad

    wcarad           = 1.25

Solvent probe radius for WCA method. This specify the radius of solvent molecules and used to define various solvent-related surfaces and volumes. This keyword is ignored when srfm is `spl2`. **For more details, see apolar keyword `srad` [here][apolar-keywords]{:target="_blank"}**

***

##### bconc

    bconc            =  0.033428

This parameter specifies the bulk solvent density in Å<sup>-3</sup>.   This can set to zero to eliminate integral contributions to nonpolar solvation energy calculation. **For more details, see apolar keyword `bconc` [here][apolar-keywords]{:target="_blank"}**.

***

##### dpos

    dpos             = 0.05

This specify the displacement in Å of the atomic positions for surface area derivative calculation using finite difference method. **For more details, see apolar keyword `dpos` [here][apolar-keywords]{:target="_blank"}**.

***

##### grid

    grid             = 0.4 0.4 0.4

This parameter specifies the quadrature grid spacing in Å for volume integral calculations. **For more details, see apolar keyword `grid` [here][apolar-keywords]{:target="_blank"}**.

***

##### APsdens

    APsdens          = 200

This specifies the number of quadrature points per Å<sup>2</sup>  to use for molecular surface or solvent accessible surface. The parameter is ignored when the `srad = 0.0` or when `srfm = spl2`. **For more details, see apolar keyword `sdens` [here][apolar-keywords]{:target="_blank"}**.

***

##### APsrfm

    APsrfm           = sacc

This parameter specifies the model used to construct the solvent-related surface and volume. **For more details, see apolar keyword `srfm` [here][apolar-keywords]{:target="_blank"}**.

***

##### APswin

    APswin           = 0.3

This specify the size of spline window in Å for  spline-based surface definitions.**For more details, see apolar keyword `swin` [here][apolar-keywords]{:target="_blank"}**.

***

##### APtemp

    APtemp           = 300 

This parameter specifies the temperature for apolar calculation. **For more details, see apolar keyword `temp` [here][apolar-keywords]{:target="_blank"}**.

***


[apolar-keywords]: http://www.poissonboltzmann.org/docs/apbs-overview/#apolar