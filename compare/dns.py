import pickle
#re_tau* = 550
filename = 'purdue/DNS_Channel_dimensional_M6.0_Reb20000_500x250x300.p'
objects = []

with (open(filename, "rb")) as openfile:
	while True:
		try:
			objects.append(pickle.load(openfile))
		except EOFError:
			break
dns = objects[0]
#print(dns['M6.0'])
print(dns['M6.0']['Data'].keys())
print(dns['M6.0']['Data']['TKE_Budget_Favre'][0])
print(dns['M6.0']['Data']['TKE_Budget_Aux'][0])
#print(dns['M6.0']['Data']['First_and_Second_Order_Mean_Favre'])
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

for i in range(N):
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
