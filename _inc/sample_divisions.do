* sample_divisions.do -- defines sample_divisions, which picks n divisions evenly spaced through the sorted list of div_id values in memory.
* Included by: calib/calib_01_alpha_cap.do, calib_02_resolution.do, calib_03_hc_k.do, calib_05_validity.do.
* Expects: div_id in memory (a division lookup); n() from $n_sample_divs (25 in _master.do).
* Notes:  we sample systematically rather than at random so the calibration sweeps see divisions from every latitude band and the sample is identical on every run without a seed.

capture program drop sample_divisions
program define sample_divisions, rclass
    syntax , n(integer)

    quietly levelsof div_id, local(all_divs)
    local ndivs : word count `all_divs'

    local sample_divs ""
    if `ndivs' <= `n' {
        local sample_divs "`all_divs'"
    }
    else {
        local prev_pos 0
        forvalues i = 1/`n' {
            local pos = round((`i' - 0.5) * `ndivs' / `n')
            if `pos' < 1 local pos 1
            if `pos' > `ndivs' local pos `ndivs'
            if `pos' != `prev_pos' {
                local d : word `pos' of `all_divs'
                local sample_divs "`sample_divs' `d'"
                local prev_pos `pos'
            }
        }
    }

    return local sample_divs "`sample_divs'"
    return local ndivs = `ndivs'
end
