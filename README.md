# BPI
Blind Power Identification (DATE 2017)
Sherief Reda (Brown University) and Adel Belouchrani (ENP)
Based on "Blind Identification of Power Sources in Processors", in IEEE/ACM Design, Automation & Test in Europe, 2017.
sherief_reda@Brown.edu and adel.belouchrani@enp.edu.dz

Sample data from HotSpot are given in the DATA folder for a 4-core configuration.

nat_trace_*:  these are samples for natural response; each line has four measurements; temperature per core; ambient temperature already substracted

steady_state.txt has 5 columns. First column is total power measurement and the remaining four are the steady-states temperatures (one per core); ambient temperature already substracted.

eval_cores is the data traces for runtime evaluation, which includes transients. It has 8 columns; the first four columns are the power per core that was given to HotSpot and the remaining four columns are the resultant temperatures. Note that the first four columns are only used to verify the accuracy. Ambient temperature already substracted
