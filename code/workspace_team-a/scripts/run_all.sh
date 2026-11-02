#!/bin/bash
# Re-run the Q0 pipeline (outputs in /workspace/out). Needs: numpy scipy astropy pandas matplotlib; astroquery (pip --user) only for fetch_catalog.py.
set -e
cd /workspace/scripts
python3 fetch_catalog.py            # Splatalogue (CDMS+JPL) dump for every spw range  (needs internet + astroquery; not run with -I)
python3 -I cont_survey.py           # peaks in every continuum/mfs image
python3 -I spectra.py               # nuclear-position spectra of all small cubes
python3 -I spectra_xb4b.py          # same for the 1.3 GB Xb4b cubes (memmap)
python3 -I blind_search.py          # first blind search at nuclei
python3 -I sysvel.py                # systemic velocity from strong lines
python3 -I moments.py; python3 -I xb4b_moments.py   # moment-0 maps
python3 -I regions.py               # objective region definition (Xb4d peaks)
python3 -I measure_lines.py         # standard line table (all regions x lines x cubes)
python3 -I make_tables.py           # status, K km/s, limits, detection matrix
python3 -I blind_all.py             # blind search in every region spectrum + catalogue candidates
python3 -I continuum.py; python3 -I contam.py       # continuum inventory + line contamination
python3 -I beam_overlap.py; python3 -I blend_fits.py; python3 -I peaks_table.py
python3 -I plot_maps.py; python3 -I wing_maps.py; python3 -I plot_spectra.py; python3 -I plot_weak.py
# ---- turn 2 (reviewer fixes) ----
python3 -I ident.py --test; python3 -I ident_blind.py > ../out/ident_blind.log
python3 -I ch3oh.py; python3 -I ch3oh_lines.py > ../out/ch3oh_lines.log; python3 -I ch3oh_matched.py; python3 -I ch3oh_ltefit.py; python3 -I ch3oh_radex.py; python3 -I line251.py
python3 -I disk.py; python3 -I disk2.py; python3 -I disk_stack.py; python3 -I regions2.py; python3 -I wclump.py
python3 -I wing_null.py; python3 -I h13cn.py; python3 -I make_table_v2.py; python3 -I test_pipeline.py
