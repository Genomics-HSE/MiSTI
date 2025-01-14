#!/bin/bash
./utils/calc_time.py ./san_sar/HGDP01029.psmc ./san_sar/HGDP00665.psmc >./san_sar/san_sar.times.txt
./utils/calc_time.py ./din_sar/DNK02.psmc ./din_sar/HGDP00665.psmc >./din_sar/din_sar.times.txt
./utils/calc_time.py ./san_din/HGDP01029.psmc ./san_din/DNK02.psmc >./san_din/san_din.times.txt
./utils/calc_time.py ./han_fre/HGDP00778.psmc ./han_fre/HGDP00521.psmc >./han_fre/han_fre.times.txt

