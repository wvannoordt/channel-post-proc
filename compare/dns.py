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

sys.exit(0)
# print(dns['M6.0']['Data'].keys())
# print(dns['M6.0']['Data']['TKE_Budget_Favre'][0])
# print(dns['M6.0']['Data']['TKE_Budget_Aux'][0])
# print(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'])
#print(len(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1]))
#print(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][1])

N = len(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1])
#print(N)
y = []
u = []
T = []
p = []

ruu  = []
rvv  = []
rww  = []
TpTp = []
upTp = []
vpTp = []
wpTp = []

for i in range(round(N/2)):
	y.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][0])
	p.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][1] * 287.15 * dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][6])
	u.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][2])
	T.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][6])
	ruu.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][10])
	rvv.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][11])
	rww.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][12])
	TpTp.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][17])
	upTp.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][18])
	vpTp.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][19])
	wpTp.append(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'][1][i][20])

outfile1 = 'purdue/cs-u.dat'
with open(outfile1, 'w') as f:
	for j in range(len(y)):
		f.write ('{} {}\n'.format(y[j], u[j]))

outfile2 = 'purdue/cs-T.dat'
with open(outfile2, 'w') as f:
	for j in range(len(y)):
		f.write ('{} {}\n'.format(y[j], T[j]))

with open('purdue/cs-ruu.dat', 'w') as f:
	for j in range(len(y)):
		f.write ('{} {}\n'.format(y[j], ruu[j]))

with open('purdue/cs-rvv.dat', 'w') as f:
	for j in range(len(y)):
		f.write ('{} {}\n'.format(y[j], rvv[j]))
		
with open('purdue/cs-rww.dat', 'w') as f:
	for j in range(len(y)):
		f.write ('{} {}\n'.format(y[j], rww[j]))

with open('purdue/cs-p.dat', 'w') as f:
	for j in range(len(y)):
		f.write ('{} {}\n'.format(y[j], p[j]))

with open('purdue/cs-TpTp.dat', 'w') as f:
	for j in range(len(y)):
		f.write ('{} {}\n'.format(y[j], TpTp[j]))

with open('purdue/cs-upTp.dat', 'w') as f:
	for j in range(len(y)):
		f.write ('{} {} {} {}\n'.format(y[j], upTp[j], vpTp[j], wpTp[j]))
