import sys
import corrfit

#c0 = sys.argv[1]     # main correlator



corrfit.all("corrfit_input_E52.txt", rdir="Post04", run="s5m90", method="gau")



