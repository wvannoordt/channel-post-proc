import sys
import pickle
#re_tau* = 550
filename = 'purdue/DNS_Channel_dimensional_M6.0_Reb20000_500x250x300.p'
objects = []

Main_File_Name = 'purdue/DNS_Channel_dimensional_M6.0_Reb20000_500x250x300.p'
f = open(Main_File_Name, 'rb')
CHANNEL_DATA = pickle.load(f,encoding='latin1')
f.close()
	
# with (open(filename, 'rb')) as openfile:
# 	while True:
# 		try:
# 			objects.append(pickle.load(openfile,encoding='latin1'))
# 		except EOFError:
# 			break
dns = CHANNEL_DATA
print(dns['M6.0']['Data'].keys())
print(dns['M6.0']['Data']['bar_A'][0])
print(dns['M6.0']['Data']['bar_AB'][0])
print(dns['M6.0']['Data']['bar_ABC'][0])
print(dns['M6.0']['Data']['bar_AB_gradC'][0])
print(dns['M6.0']['Data']['bar_A_gradB'][0])
print(dns['M6.0']['Data']['bar_rhoAppBpp'][0])
print(dns['M6.0']['Data']['bar_dtau_ij_dxl'][0]) # i first, then j, then l
print(dns['M6.0']['Data']['bar_dtau_ij_dxl'][1].shape) # i first, then j, then l
# (1, 250, 9, 3)
#  ^  ^    ^  ^ l
#  |  |    | i,j -> (i,j) -> (i-1)*3 + j
#  |  |
#  | redundant
print(dns['M6.0']['Data']['bar_tau_ij'][0])
print(dns['M6.0']['Data']['bar_ui_dtau_jl_dxm'][0]) # i first, then j, then l, then m
print(dns['M6.0']['Data']['bar_ui_dtau_jl_dxm'][1].shape) # i first, then j, then l, then m

# (1, 250, 9, 9)
#  ^       ^  ^ l,m -> (l,m) -> (l-1)*3 + m
#  |       | i,j -> (i,j) -> (i-1)*3 + j
#  |  ^
#  | redundant

A00 = []
A01 = []
A02 = []
A10 = []
A11 = []
A12 = []
A20 = []
A21 = []
A22 = []
