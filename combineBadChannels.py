#! /usr/bin/env python
from os import path
import sys

# import bad_channel_lhc15o_1
# import bad_channel_lhc15o_2
# import bad_channel_lhc15o_3
# import bad_channel_lhc15o_5

import bad_channel_lhc17q_all1
import bad_channel_lhc17q_all2
#import bad_channel_lhc15o_3_new
#import bad_channel_lhc15o_4_new

# in2not3 = []
# bad_all1 = sorted(bad_channel_lhc15o_1.bad_all)
# bad_all2 = sorted(bad_channel_lhc15o_2.bad_all)
# bad_all3 = sorted(bad_channel_lhc15o_3.bad_all)
# bad_all5 = sorted(bad_channel_lhc15o_5.bad_all)

bad_all1_new = sorted(bad_channel_lhc17q_all1.bad_all)
bad_all2_new = sorted(bad_channel_lhc17q_all2.bad_all)
#bad_all3_new = sorted(bad_channel_lhc15o_3_new.bad_all)
#bad_all4_new = sorted(bad_channel_lhc15o_4_new.bad_all)

bad_Allchn = [bad_all1_new,bad_all2_new]#, bad_all3_new, bad_all4_new]
bad_all1234 = []

#####################################################################################################################

dead1 = sorted(bad_channel_lhc17q_all1.dead_emcal + bad_channel_lhc17q_all1.dead_dcal)
dead2 = sorted(bad_channel_lhc17q_all2.dead_emcal + bad_channel_lhc17q_all2.dead_dcal)
#dead3 = sorted(bad_channel_lhc17q_all3.dead_emcal + bad_channel_lhc17q_all3.dead_dcal)
#dead4 = sorted(bad_channel_lhc17q_all4.dead_emcal + bad_channel_lhc17q_all4.dead_dcal)

dead_chn = [dead1_new,dead2_new]#, dead3_new, dead4_new]
dead1234 = []


bad1 = sorted(bad_channel_lhc17q_all1.bad_emcal + bad_channel_lhc17q_all1.bad_dcal)
bad2 = sorted(bad_channel_lhc17q_all2.bad_emcal + bad_channel_lhc17q_all2.bad_dcal)
#bad3 = sorted(bad_channel_lhc17q_all3.bad_emcal + bad_channel_lhc17q_all3.bad_dcal)
#bad4 = sorted(bad_channel_lhc17q_all4.bad_emcal + bad_channel_lhc17q_all4.bad_dcal)

bad_chn = [bad1_new,bad2_new]#, bad3_new, bad4_new]
bad1234 = []


warm1 = sorted(bad_channel_lhc17q_all1.warm_emcal + bad_channel_lhc17q_all1.warm_dcal)
warm2 = sorted(bad_channel_lhc17q_all2.warm_emcal + bad_channel_lhc17q_all2.warm_dcal)
#warm3 = sorted(bad_channel_lhc17q_all3.warm_emcal + bad_channel_lhc17q_all3.warm_dcal)
#warm4 = sorted(bad_channel_lhc17q_all4.warm_emcal + bad_channel_lhc17q_all4.warm_dcal)

warm_chn = [warm1_new,warm2_new]#, warm3_new, warm4_new]
warm1234 = []



outliers = []
for i in bad_all1_new:
    
