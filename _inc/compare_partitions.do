* compare_partitions.do -- defines compare_partitions, which returns agreement statistics (MI, NMI, ARI, AMI, percent identical) between two partition variables in memory.
* Included by: assess/assess_02_agreement.do, assess_03_temporal.do, assess_04_sensitivity.do, assess_06_czones.do, _inc/compare_to_baseline.do.
* Expects: the a() and b() partition variables, one row per unit; $py and $program for the sidecar.
* Reads:  the sidecar's output CSV (tempfile)
* Writes: a tempfile CSV of a() b() for shell "$py" _py/compare_partitions.py
* Notes:  MI, NMI (arithmetic normalisation) and ARI are computed in Stata from the contingency table; AMI and NMI_geom come from scikit-learn through the sidecar. r(ARI_py) is the sidecar's ARI, returned for cross-checking against r(ARI).

capture program drop compare_partitions
program define compare_partitions, rclass
    syntax , a(name) b(name)

    quietly count
    local N = r(N)
    return scalar N = `N'

    local AMI       = .
    local AMI_ari   = .
    local NMI_geom  = .
    capture confirm file "$program/_py/compare_partitions.py"
    if !_rc {
        tempfile pbase
        local incsv  "`pbase'_in.csv"
        local outcsv "`pbase'_out.csv"
        export delimited `a' `b' using "`incsv'", replace nolabel
        shell "$py" "$program/_py/compare_partitions.py" "`incsv'" "`outcsv'"
        capture confirm file "`outcsv'"
        if !_rc {
            tempname pf
            frame create `pf'
            frame `pf' {
                quietly import delimited using "`outcsv'", varnames(1) clear
                local AMI      = ami[1]
                local AMI_ari  = ari[1]
                local NMI_geom = nmi_geom[1]
            }
            frame drop `pf'
        }
        else {
            display as error "compare_partitions: AMI sidecar produced no output; AMI set to ."
        }
    }
    return scalar AMI      = `AMI'
    return scalar NMI_geom = `NMI_geom'
    return scalar ARI_py   = `AMI_ari'

    quietly count if `a' == `b'
    return scalar pct_stable = r(N) / `N' * 100

    preserve
    contract `a' `b'
    rename _freq n_ij

    bysort `a': egen long a_i = total(n_ij)
    bysort `b': egen long b_j = total(n_ij)

    gen double mi_term = (n_ij / `N') * log(`N' * n_ij / (a_i * b_j)) if n_ij > 0
    quietly summarize mi_term
    local MI = r(sum)

    gen double c_n2 = n_ij * (n_ij - 1) / 2
    quietly summarize c_n2
    local Index = r(sum)

    tempfile contingency
    save `contingency'

    keep `a' a_i
    duplicates drop
    gen double h_term = -(a_i / `N') * log(a_i / `N') if a_i > 0
    quietly summarize h_term
    local HX = r(sum)
    gen double c_a2 = a_i * (a_i - 1) / 2
    quietly summarize c_a2
    local SUM_A = r(sum)

    use `contingency', clear

    keep `b' b_j
    duplicates drop
    gen double h_term = -(b_j / `N') * log(b_j / `N') if b_j > 0
    quietly summarize h_term
    local HY = r(sum)
    gen double c_b2 = b_j * (b_j - 1) / 2
    quietly summarize c_b2
    local SUM_B = r(sum)

    local NMI = cond(`HX' + `HY' > 0, 2 * `MI' / (`HX' + `HY'), 0)

    local CN2      = `N' * (`N' - 1) / 2
    local Expected = `SUM_A' * `SUM_B' / `CN2'
    local Max      = 0.5 * (`SUM_A' + `SUM_B')
    local ARI      = cond(`Max' - `Expected' != 0, ///
        (`Index' - `Expected') / (`Max' - `Expected'), 0)

    return scalar MI  = `MI'
    return scalar NMI = `NMI'
    return scalar ARI = `ARI'
    restore
end
