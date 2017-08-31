#! /usr/bin/env python
from os import path
import sys

# import bad_channel_lhc15o_1
# import bad_channel_lhc15o_2
# import bad_channel_lhc15o_3
# import bad_channel_lhc15o_5

import bad_channel_lhc15o_1_new
import bad_channel_lhc15o_2_new
import bad_channel_lhc15o_3_new
import bad_channel_lhc15o_4_new

# in2not3 = []
# bad_all1 = sorted(bad_channel_lhc15o_1.bad_all)
# bad_all2 = sorted(bad_channel_lhc15o_2.bad_all)
# bad_all3 = sorted(bad_channel_lhc15o_3.bad_all)
# bad_all5 = sorted(bad_channel_lhc15o_5.bad_all)

bad_all1_new = sorted(bad_channel_lhc15o_1_new.bad_all)
bad_all2_new = sorted(bad_channel_lhc15o_2_new.bad_all)
bad_all3_new = sorted(bad_channel_lhc15o_3_new.bad_all)
bad_all4_new = sorted(bad_channel_lhc15o_4_new.bad_all)

bad_chn = [bad_all1_new,bad_all2_new, bad_all3_new, bad_all4_new]
bad_all1234 = []
        
for a in bad_chn:
    for n in a:
        if not (n in bad_all1234):
            bad_all1234.append(n)
            
outliers = [1824, 1826, 1827, 1829, 1830, 1831, 1833, 1834, 1835, 1872, 1874, 1878, 1880, 1882, 1884, 1886, 3176, 7507, 12480, 12481, 12482, 12483, 12484, 12485, 12486, 12488, 12489, 12490, 12491, 12492, 12493, 12494, 12495, 12528, 12530, 12531, 12532, 12533, 12534, 12535, 12536, 12537, 12538, 12539, 12540, 12541, 12542, 12543, 14240, 14241, 14242, 14243, 14244, 14245, 14246, 14248, 14250, 14251, 14252, 14253, 14288, 14289, 14290, 14291, 14292, 14293, 14294, 14295, 14296, 14297, 14298, 14299, 14300, 14301, 14384, 14386, 14388, 14390, 14392, 14394]

for n in outliers:
    if n not in bad_all1234:
        print n, 'is not in bad_all1234'
