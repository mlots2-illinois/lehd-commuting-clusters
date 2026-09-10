* build_div_dissim.do -- defines build_div_dissim, which builds the within-division pair table (distance, P_ij and the dissimilarity D) from a division lookup and a pooled flow file.
* Included by: calib/calib_01_alpha_cap.do, calib/calib_02_resolution.do, calib/calib_03_hc_k.do; cluster/cluster_01_hc.do includes it but runs _py/build_dissim.py instead.
* Expects: the divlookup() and pflow() files named by the caller; outdir() also needs alpha(), cap_km() and divfmt().
* Reads:  divlookup() (div_id geoid lat lon); pflow() (geo_i geo_j P_ij)
* Writes: <outdir>/dissim_div<NNNNN>.csv (one per division; stale files erased first) and/or the pairsout() .dta of all pairs
* Notes:  D = alpha * (1 - min(P_ij, 1)) + (1 - alpha) * min(d_km / cap_km, 1), the same formula as _py/build_dissim.py. A pair with no flow row gets P_ij = 0. The calibration sweeps use this version on sampled divisions; the full run uses the Python twin.

capture program drop build_div_dissim
program define build_div_dissim
    syntax , divlookup(string) pflow(string)                          ///
        [alpha(string) cap_km(string) outdir(string)             ///
        divfmt(string) restrict_divs(string) keepvars           ///
        pairsout(string) countsframe(name)]

    if "`outdir'" == "" & "`pairsout'" == "" {
        display as error "build_div_dissim: specify outdir() and/or pairsout()."
        exit 198
    }
    if "`outdir'" != "" {
        if "`alpha'" == "" | "`cap_km'" == "" | "`divfmt'" == "" {
            display as error ///
                "build_div_dissim: outdir() requires alpha(), cap_km(), and divfmt()."
            exit 198
        }
    }
    local do_D = ("`alpha'" != "" & "`cap_km'" != "")

    if "`outdir'" != "" {
        capture mkdir "`outdir'"

        local _stale : dir "`outdir'" files "dissim_div*.csv"
        foreach f of local _stale {
            capture erase "`outdir'/`f'"
        }
    }

    capture frame drop _bdd_pflow
    frame create _bdd_pflow
    frame _bdd_pflow: use geo_i geo_j P_ij using "`pflow'", clear

    capture frame drop _bdd_units
    frame create _bdd_units
    frame _bdd_units: use div_id geoid lat lon using "`divlookup'", clear
    if "`restrict_divs'" != "" {
        frame _bdd_units {
            gen byte _keep = 0
            foreach d of local restrict_divs {
                replace _keep = 1 if div_id == `d'
            }
            keep if _keep == 1
            drop _keep
        }
    }

    frame _bdd_units: quietly levelsof div_id, local(divs)

    if "`countsframe'" != "" {
        capture frame drop `countsframe'
        frame create `countsframe' long div_id long n_units
    }

    local pair_pieces ""

    foreach d of local divs {

        capture frame drop _bdd_d
        frame _bdd_units: frame put geoid lat lon if div_id == `d', into(_bdd_d)

        frame _bdd_d: local N = _N
        if `N' < 2 {
            display as text "  build_div_dissim div `d': only `N' unit(s), skipping"
            continue
        }
        if "`countsframe'" != "" {
            frame post `countsframe' (`d') (`N')
        }

        display as text "  build_div_dissim div `d': `N' units, building ~" ///
            %12.0fc (`N' * (`N' - 1) / 2) " pairs ..."

        frame _bdd_d {

            sort geoid

            tempfile div_units
            save `div_units'

            rename (geoid lat lon) (geo_j lat_j lon_j)
            tempfile jcopy
            save `jcopy'

            use `div_units', clear
            rename (geoid lat lon) (geo_i lat_i lon_i)
            cross using `jcopy'
            keep if geo_i < geo_j

            quietly {
                local torad = _pi / 180
                gen double a = sin((lat_j - lat_i)*`torad'/2)^2 ///
                    + cos(lat_i*`torad')*cos(lat_j*`torad') ///
                    *sin((lon_j - lon_i)*`torad'/2)^2
                gen double d_km = 2 * 6371 * atan2(sqrt(a), sqrt(1 - a))
                drop a

                frlink m:1 geo_i geo_j, frame(_bdd_pflow) generate(_bdd_lnk)
                frget P_ij, from(_bdd_lnk)
                drop _bdd_lnk
                replace P_ij = 0 if mi(P_ij)
            }

            if `do_D' {
                gen double P_clip = min(P_ij, 1)
                gen double D = `alpha' * (1 - P_clip) ///
                    + (1 - `alpha') * min(d_km / `cap_km', 1)
                drop P_clip
            }

            sort geo_i geo_j

            if "`outdir'" != "" {
                local fname = string(`d', "`divfmt'")
                if "`keepvars'" != "" {
                    quietly export delimited geo_i geo_j d_km P_ij D ///
                        using "`outdir'/dissim_div`fname'.csv", replace
                }
                else {
                    quietly export delimited geo_i geo_j D ///
                        using "`outdir'/dissim_div`fname'.csv", replace
                }
            }

            if "`pairsout'" != "" {
                gen long div_id = `d'
                if `do_D' {
                    keep div_id geo_i geo_j d_km P_ij D
                    order div_id geo_i geo_j d_km P_ij D
                }
                else {
                    keep div_id geo_i geo_j d_km P_ij
                    order div_id geo_i geo_j d_km P_ij
                }
                tempfile _bdp_`d'
                quietly save "`_bdp_`d''"
                local pair_pieces `"`pair_pieces' "`_bdp_`d''""'
            }

            display as text "  build_div_dissim div `d': `N' units, " %9.0fc _N " pairs"
        }
    }

    if "`pairsout'" != "" {
        clear
        if `"`pair_pieces'"' != "" {
            append using `pair_pieces'
        }
        sort div_id geo_i geo_j
        compress
        save "`pairsout'", replace
        display as text "  build_div_dissim: pairs saved -> `pairsout' (" _N " rows)"
    }

    capture frame drop _bdd_d
    capture frame drop _bdd_units
    capture frame drop _bdd_pflow
end
