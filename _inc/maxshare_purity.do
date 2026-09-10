* maxshare_purity.do -- defines maxshare_purity, which collapses to one row per reference unit with the share of its members that sit in its largest cluster.
* Included by: assess/assess_05_cbsa.do, assess/assess_06_czones.do.
* Expects: one row per unit in memory with the refvar() and clustervar() variables.
* Notes:  purity = max_c n(ref, c) / n(ref). altgen() gives 1 / (number of clusters in the reference unit), the value purity would take if the members were spread evenly. The data in memory are replaced.

capture program drop maxshare_purity
program define maxshare_purity
    version 18
    syntax , REFvar(varname) CLUSTERvar(varname) GEN(name) [ALTGen(name)]

    confirm new variable `gen'
    if "`altgen'" != "" confirm new variable `altgen'

    bysort `refvar' `clustervar': gen long _msp_nuc = _N
    by     `refvar' `clustervar': keep if _n == 1

    bysort `refvar': egen long _msp_nc  = total(_msp_nuc)
    bysort `refvar': egen long _msp_mx  = max(_msp_nuc)
    bysort `refvar': gen  long _msp_ncl = _N
    by     `refvar': keep if _n == 1

    gen double `gen' = _msp_mx / _msp_nc
    if "`altgen'" != "" {
        gen double `altgen' = 1 / _msp_ncl
    }

    drop _msp_nuc _msp_nc _msp_mx _msp_ncl
end
