import sys
import pickle
#re_tau* = 550
def main():
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
	print(dns['M6.0']['Data']['bar_tau_ij_dul_dxm'][0]) # i first, then j, then l
	print(dns['M6.0']['Data']['bar_tau_ij_dul_dxm'][1].shape) # i first, then j, then l
	print(dns['M6.0']['Data']['bar_gradA'][0]) # i first, then j, then l

	# (1, 250, 9, 3)
	#  ^  ^    ^  ^ l
	#  |  |    | i,j -> (i,j) -> (i-1)*3 + j
	#  |  |
	#  | redundant
	print(dns['M6.0']['Data']['bar_tau_ij'][0])
	print(dns['M6.0']['Data']['bar_ui_dtau_jl_dxm'][0]) # i first, then j, then l, then m
	print(dns['M6.0']['Data']['bar_ui_dtau_jl_dxm'][1].shape) # i first, then j, then l, then m
	print(dns['M6.0']['Data']['y'][0]) # i first, then j, then l, then m

	# (1, 250, 9, 9)
	#  ^       ^  ^ l,m -> (l,m) -> (l-1)*3 + m
	#  |       | i,j -> (i,j) -> (i-1)*3 + j
	#  |  ^
	#  | redundant
	N = dns['M6.0']['Data']['bar_ui_dtau_jl_dxm'][1].shape[0]

	y   = []

	A00 = []
	A01 = []
	A02 = []
	A10 = []
	A11 = []
	A12 = []
	A20 = []
	A21 = []
	A22 = []

	D0  = []
	D1  = []
	D2  = []

	B00 = []
	B01 = []
	B02 = []

	C10 = []

	u_bar = []
	u_tld = []
	v_bar = []
	v_tld = []
	w_bar = []
	w_tld = []
	T_bar = []
	T_tld = []
	
	u_pp  = []
	v_pp  = []
	w_pp  = []
	
	uT_pp  = []
	vT_pp  = []
	wT_pp  = []


	tau_du = dns['M6.0']['Data']['bar_tau_ij_dul_dxm'][1]
	tau    = dns['M6.0']['Data']['bar_tau_ij'][1]
	grads  = dns['M6.0']['Data']['bar_gradA'][1]
	vals   = dns['M6.0']['Data']['bar_A'][1]
	vals2  = dns['M6.0']['Data']['bar_AB'][1]
	dtau   = dns['M6.0']['Data']['bar_dtau_ij_dxl'][1]
	yy     = dns['M6.0']['Data']['y'][1]
	trips  = dns['M6.0']['Data']['bar_ABC'][1]
	pps    = dns['M6.0']['Data']['bar_rhoAppBpp'][1]

	# print(grads)
	print(N, round(N/2))
	for i in range(round(N/2)):
		y.append(yy[i])
		u_bar.append(vals[i][1])
		v_bar.append(vals[i][2])
		w_bar.append(vals[i][3])
		T_bar.append(vals[i][5])
		u_tld.append(vals2[i][0]/vals[i][0])
		v_tld.append(vals2[i][1]/vals[i][0])
		w_tld.append(vals2[i][2]/vals[i][0])
		T_tld.append(vals2[i][3]/vals[i][0])
		A00.append(tau_du[i][0][0][0][0] - tau[i][0][0]*grads[i][1][0])
		A01.append(tau_du[i][0][1][0][1] - tau[i][0][1]*grads[i][1][1])
		A02.append(tau_du[i][0][2][0][2] - tau[i][0][2]*grads[i][1][2])
		A10.append(tau_du[i][1][0][1][0] - tau[i][1][0]*grads[i][2][0])
		A11.append(tau_du[i][1][1][1][1] - tau[i][1][1]*grads[i][2][1])
		A12.append(tau_du[i][1][2][1][2] - tau[i][1][2]*grads[i][2][2])
		A20.append(tau_du[i][2][0][2][0] - tau[i][2][0]*grads[i][3][0])
		A21.append(tau_du[i][2][1][2][1] - tau[i][2][1]*grads[i][3][1])
		A22.append(tau_du[i][2][2][2][2] - tau[i][2][2]*grads[i][3][2])
		D0.append (vals[i][4]*grads[i][1][0] - 0)
		D1.append (vals[i][4]*grads[i][2][1] - 0)
		D2.append (vals[i][4]*grads[i][3][2] - 0)
		B00.append((u_tld[i] - u_bar[i])*dtau[i][0][0][0])
		B01.append((u_tld[i] - u_bar[i])*dtau[i][0][1][1])
		B02.append((u_tld[i] - u_bar[i])*dtau[i][0][2][2])
		C10.append((u_bar[i] - u_tld[i])*grads[i][4][1])
		u_pp.append((trips[i][0]/vals[i][0]) - (vals2[i][0]/vals[i][0])*(vals2[i][0]/vals[i][0]))
		v_pp.append((trips[i][3]/vals[i][0]) - (vals2[i][1]/vals[i][0])*(vals2[i][1]/vals[i][0]))
		w_pp.append((trips[i][5]/vals[i][0]) - (vals2[i][2]/vals[i][0])*(vals2[i][2]/vals[i][0]))
		uT_pp.append(pps[i][6]/vals[i][0])
		vT_pp.append(pps[i][7]/vals[i][0])
		wT_pp.append(pps[i][8]/vals[i][0])

	write_csv(y, A00,   'purdue/cs-a00.csv')
	write_csv(y, A01,   'purdue/cs-a01.csv')
	write_csv(y, A02,   'purdue/cs-a02.csv')
	write_csv(y, A10,   'purdue/cs-a10.csv')
	write_csv(y, A11,   'purdue/cs-a11.csv')
	write_csv(y, A12,   'purdue/cs-a12.csv')
	write_csv(y, A20,   'purdue/cs-a20.csv')
	write_csv(y, A21,   'purdue/cs-a21.csv')
	write_csv(y, A22,   'purdue/cs-a22.csv')
	write_csv(y, D0,    'purdue/cs-d0.csv')
	write_csv(y, D1,    'purdue/cs-d1.csv')
	write_csv(y, D2,    'purdue/cs-d2.csv')
	write_csv(y, B00,   'purdue/cs-b00.csv')
	write_csv(y, B01,   'purdue/cs-b01.csv')
	write_csv(y, B02,   'purdue/cs-b02.csv')
	write_csv(y, C10,   'purdue/cs-c10.csv')
	write_csv(y, T_tld, 'purdue/cs-T.csv')
	write_csv(y, u_tld, 'purdue/cs-u.csv')
	write_csv(y, u_pp, 'purdue/cs-upp.csv')
	write_csv(y, v_pp, 'purdue/cs-vpp.csv')
	write_csv(y, w_pp, 'purdue/cs-wpp.csv')
	write_csv(y, uT_pp, 'purdue/cs-uTpp.csv')
	write_csv(y, vT_pp, 'purdue/cs-vTpp.csv')
	write_csv(y, wT_pp, 'purdue/cs-wTpp.csv')

def write_csv(y, f, name):
	print('output {}'.format(name))
	with open(name, 'w') as fh:
		for i in range(len(y)):
			fh.write('{},{}\n'.format(y[i], f[i]))

if __name__ == "__main__":
	main()
	sys.exit(0)