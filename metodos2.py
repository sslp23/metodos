from sympy import *
x, y, z, t = symbols('x y z t')

arquivo = open('saida.txt', 'w')


def solver(func, yn, tn):
	return func.subs({y:yn, t:tn})

def euler(y0, t0, h, n, func):
	print("y( 0.0 ) = ", y0, file=arquivo)
	print("h = ", h, file=arquivo)
	print("0 ", y0, file=arquivo)

	y = y0
	t = t0
	for i in range(0, n, 1):
		k = solver(func, y, t)
		y = y + h*k
		t+=h
		print(i+1, ' ', y, file=arquivo)

	return

def euler_list(y0, t0, h, n, func): #funcao q salva uma lista de valores achados por euler p metodos de passo multiplo
	y = list(range(n))
	y[0] = y0
	t = t0
	for i in range(0, n-1, 1):
		k = solver(func, y[i], t)
		y[i+1] = y[i] + h*k
		t+=h

	return y

def euler_inverso(y0, t0, h, n, func):
	print("y( 0.0 ) = ", y0, file=arquivo)
	print("h = ", h, file=arquivo)
	print("0 ", y0, file=arquivo)

	yn = y0
	tn = t0
	for i in range(0, n, 1):
		tn+=h
		k = yn + solver(func, y, tn)*h
		sol = solve(Eq(k, y), y)
		yn = sol[0]
		print(i+1, ' ', yn, file=arquivo)
	
	return

def euler_inverso_prev(y0, t0, h, n, func):
	print("y( 0.0 ) = ", y0, file=arquivo)
	print("h = ", h, file=arquivo)
	print("0 ", y0, file=arquivo)

	yn = y0
	tn = t0
	k1 = 0
	for i in range(0, n, 1):
		k1 = yn + solver(func, yn, tn)*h #previsao com euler
		tn+=h
		yn = yn + solver(func, k1, tn)*h
		print(i+1, ' ', yn, file=arquivo)
	
	return

def euler_inverso_list(y0, t0, h, n, func):
	y = list(range(n))
	y[0] = y0
	tn = t0
	k1 = 0
	for i in range(0, n-1, 1):
		k1 = y[i] + solver(func, y[i], tn)*h #previsao com euler
		tn+=h
		y[i+1] = y[i] + solver(func, k1, tn)*h
	
	return y

def euler_aprimorado(y0, t0, h, n, func):
	print("y( 0.0 ) = ", y0, file=arquivo)
	print("h = ", h, file=arquivo)
	print("0 ", y0, file=arquivo)

	y = y0;
	t = t0;
	for i in range(0, n, 1):
		k1 = solver(func, y, t)
		k2 = solver(func, y+(h*k1), t+h)
		y = y+((h/2)*(k1+k2))
		t += h	
		print(i+1, ' ', y, file=arquivo)

	return

def euler_aprimorado_list(y0, t0, h, n, func):
	y = list(range(n))
	y[0] = y0
	t = t0
	for i in range(0, n-1, 1):
		k1 = solver(func, y[i], t)
		k2 = solver(func, y[i]+(h*k1), t+h)
		y[i+1] = y[i]+((h/2)*(k1+k2))
		t += h	

	return y

def runge_kutta(y0, t0, h, n, func):
	print("y( 0.0 ) = ", y0, file=arquivo)
	print("h = ", h, file=arquivo)
	print("0 ", y0, file=arquivo)

	y = y0;
	t = t0;
	for i in range(0, n, 1):
		k1 = solver(func, y, t)
		k2 = solver(func, y+(h*k1*0.5), t+(h*0.5))
		k3 = solver(func, y+(h*k2*0.5), t+(h*0.5))
		k4 = solver(func, y+(h*k3), t+h)
		y = y+((h/6)*(k1+(2*k2)+(2*k3)+k4))
		t += h	
		print(i+1, ' ', y, file=arquivo)

	return

def runge_kutta_list(y0, t0, h, n, func):
	y = list(range(n))
	y[0] = y0;
	t = t0;
	for i in range(0, n-1, 1):
		k1 = solver(func, y[i], t)
		k2 = solver(func, y[i]+(h*k1*0.5), t+(h*0.5))
		k3 = solver(func, y[i]+(h*k2*0.5), t+(h*0.5))
		k4 = solver(func, y[i]+(h*k3), t+h)
		y[i+1] = y[i]+((h/6)*(k1+(2*k2)+(2*k3)+k4))
		t += h	

	return y

def adam_bashfort(y, t0, h, n, func, ordem):
	print("y( 0.0 ) = ", y[0], file=arquivo)
	print("h = ", h, file=arquivo)
	
	#vetor do valor do passo sendo criado:
	t = list(range(n+2))
	t[0] = t0	
	for i in range(1, n+2, 1):
		t[i] = t[i-1] + h;

	#vetor de f
	f = list(range(ordem))
	for i in range(0, ordem, 1):
		f[i] = solver(func, y[i], t[i])

	for i in range(0, ordem, 1):
		print(i, ' ', y[i], file=arquivo)

	if ordem == 2:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			y.append(y[x]+((3/2)*h*f[x])-((1/2)*h*f[x-1]))
			f.append(solver(func, y[x+1], t[x+2]))
			print(i, ' ', y[i], file=arquivo)
	
	if ordem == 3:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			y.append(y[x]+((23/12)*h*f[x])-((4/3)*h*f[x-1])+((5/12)*h*f[x-2]))
			f.append(solver(func, y[x+1], t[x+2]))
			print(i, ' ', y[i], file=arquivo)

	if ordem == 4:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			y.append(y[x]+((55/24)*h*f[x])-((59/24)*h*f[x-1])+((37/24)*h*f[x-2])-((3/8)*h*f[x-3]))
			f.append(solver(func, y[x+1], t[x+2]))
			print(i, ' ', y[i], file=arquivo)

	if ordem == 5:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			y.append(y[x]+(h*((1901/720)*f[x]-(1387/360)*f[x-1]+(109/30)*f[x-2]-(637/360)*f[x-3]+(251/720)*f[x-4])))
			f.append(solver(func, y[x+1], t[x+2]))
			print(i, ' ', y[i], file=arquivo)
	
	if ordem == 6:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			y.append(y[x]+(4277/1440)*h*f[x]-(2641/480)*h*f[x-1]+(4991/720)*h*f[x-2]-(3649/720)*h*f[x-3]+(959/480)*h*f[x-4]-(95/288)*h*f[x-5])
			f.append(solver(func, y[x+1], t[x+2]))
			print(i, ' ', y[i], file=arquivo)

	if ordem == 7:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			y.append(y[x]+(((198721/60480)*h*f[x]-(18637/2520)*h*f[x-1]+(235183/20160)*h*f[x-2]-(10754/945)*h*f[x-3]+(135713/20160)*h*f[x-4]-(5603/2520)*h*f[x-5]+(19087/60480)*h*f[x-6])))
			f.append(solver(func, y[x+1], t[x+2]))
			print(i, ' ', y[i], file=arquivo)

	if ordem == 8:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			y.append(y[x]+((((16083/4480)*h*f[x])-((1152169/120960)*h*f[x-1])+((242653/13440)*h*f[x-2])-((296053/13440)*h*f[x-3])+((2102243/120960)*h*f[x-4])-(115747/13440)*h*f[x-5]+(32863/13440)*h*f[x-6]-(5257/17280)*h*f[x-7])))
			f.append(solver(func, y[x+1], t[x+2]))
			print(i, ' ', y[i], file=arquivo)
					
	return 

#ordem = nr de pontos

def adam_moulton(y, t0, h, n, func, ordem):
	print("y( 0.0 ) = ", y[0], file=arquivo)
	print("h = ", h, file=arquivo)
	
	#vetor do valor do passo sendo criado:
	t = list(range(n+2))
	t[0] = t0	
	for i in range(1, n+2, 1):
		t[i] = t[i-1] + h;

	#vetor de f
	f = list(range(ordem))
	for i in range(0, ordem, 1):
		f[i] = solver(func, y[i], t[i])

	for i in range(0, ordem, 1):
		print(i, ' ', y[i], file=arquivo)

	if ordem == 2: 
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			yn1 = y[x]+((3/2)*h*f[x])-((1/2)*h*f[x-1])
			fn1 = solver(func, yn1, t[i]) #prevendo com adam bashfort
			y.append(y[x]+((1/2)*h*fn1)+((1/2)*h*f[x]))
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)

	if ordem == 3:
		for i in range(ordem, n+1, 1):
			yn1 = y[x]+((23/12)*h*f[x])-((4/3)*h*f[x-1])+((5/12)*h*f[x-2])
			fn1 = solver(func, yn1, t[i])#prever com adam bashfort
			x = len(y)-1
			y.append(y[x]+((5/12)*h*fn1)+((2/3)*h*f[x])-((1/12)*h*f[x-1]))
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)
	 
	if ordem == 4:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			yn1 = y[x]+((55/24)*h*f[x])-((59/24)*h*f[x-1])+((37/24)*h*f[x-2])-((3/8)*h*f[x-3])
			fn1 = solver(func, yn1, t[i])#prever com adam bashfort
			y.append(y[x]+((3/8)*h*fn1)+((19/24)*h*f[x])-((5/24)*h*f[x-1])+((1/24)*h*f[x-2]))
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)

	if ordem == 5: 
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			yn1 = y[x]+(h*((1901/720)*f[x]-(1387/360)*f[x-1]+(109/30)*f[x-2]-(637/360)*f[x-3]+(251/720)*f[x-4]))
			fn1 = solver(func, yn1, t[i])#prever com adam bashfort
			y.append(y[x]+((251/720)*h*fn1)+((323/360)*h*f[x])-((11/30)*h*f[x-1])+((53/360)*h*f[x-2])-((19/720)*h*f[x-3]))
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)

	if ordem == 6:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			yn1 = y[x]+(4277/1440)*h*f[x]-(2641/480)*h*f[x-1]+(4991/720)*h*f[x-2]-(3649/720)*h*f[x-3]+(959/480)*h*f[x-4]-(95/288)*h*f[x-5]
			fn1 = solver(func, yn1, t[i])#prever com adam bashfort
			y.append(y[x]+((95/288)*h*fn1+(1427/1440)*h*f[x]-(133/240)*h*f[x-1]+(241/720)*h*f[x-2]-(173/1440)*h*f[x-3]+(3/160)*h*f[x-4]))
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)
	
	if ordem == 7: 
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			yn1 = y[x]+(((198721/60480)*h*f[x]-(18637/2520)*h*f[x-1]+(235183/20160)*h*f[x-2]-(10754/945)*h*f[x-3]+(135713/20160)*h*f[x-4]-(5603/2520)*h*f[x-5]+(19087/60480)*h*f[x-6]))
			fn1 = solver(func, yn1, t[i])#prever com adam bashfort
			y.append(y[x]+(19087/60480)*h*fn1+(2713/2520)*h*f[x]-(15487/20160)*h*f[x-1]+(586/945)*h*f[x-2]-(6737/20160)*h*f[x-3]+(263/2520)*h*f[x-4]-(863/60480)*h*f[x-5])
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)

	if ordem == 8: #fazer
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			yn1 = y[x]+((((16083/4480)*h*f[x])-((1152169/120960)*h*f[x-1])+((242653/13440)*h*f[x-2])-((296053/13440)*h*f[x-3])+((2102243/120960)*h*f[x-4])-(115747/13440)*h*f[x-5]+(32863/13440)*h*f[x-6]-(5257/17280)*h*f[x-7]))
			fn1 = solver(func, yn1, t[i]) #prevendo com adam bashfort
			y.append(y[x]+(5257/17280)*h*fn1+(139849/120960)*h*f[x]-(4511/4480)*h*f[x-1]+(123133/120960)*h*f[x-2]-(88547/120960)*h*f[x-3]+(1537/4480)*h*f[x-4]-(11351/120960)*h*f[x-5]+(275/24192)*h*f[x-6])
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)

					
	return

def formula_inversa(y, t0, h, n, func, ordem):
	print("y( 0.0 ) = ", y[0], file=arquivo)
	print("h = ", h, file=arquivo)
	
	#vetor do valor do passo sendo criado:
	t = list(range(n+2))
	t[0] = t0	
	for i in range(1, n+2, 1):
		t[i] = t[i-1] + h;

	#vetor de f
	f = list(range(ordem))
	for i in range(0, ordem, 1):
		f[i] = solver(func, y[i], t[i])

	for i in range(0, ordem, 1):
		print(i, ' ', y[i], file=arquivo)

	if ordem == 2: 
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			yn1 = y[x]+((3/2)*h*f[x])-((1/2)*h*f[x-1])
			fn1 = solver(func, yn1, t[i]) #prevendo com adam bashfort
			y.append((4/3)*y[x]-(1/3)*y[x-1]+(2/3)*h*fn1)
			f.append(solver(func, y[x+1], t[i])) # fazendo o bashfort
			print(i, ' ', y[i], file=arquivo)

	if ordem == 3:
		for i in range(ordem, n+1, 1):
			yn1 = y[x]+((23/12)*h*f[x])-((4/3)*h*f[x-1])+((5/12)*h*f[x-2])
			fn1 = solver(func, yn1, t[i])#prever com adam bashfort
			x = len(y)-1
			y.append((18/11)*y[x]-(9/11)*y[x-1]+(2/11)*y[x-2]+(6/11)*h*fn1)
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)
	 
	if ordem == 4:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			yn1 = y[x]+((55/24)*h*f[x])-((59/24)*h*f[x-1])+((37/24)*h*f[x-2])-((3/8)*h*f[x-3])
			fn1 = solver(func, yn1, t[i])#prever com adam bashfort
			y.append((48/25)*y[x]-(36/25)*y[x-1]+(16/25)*y[x-2]-(3/25)*y[x-3]+(12/25)*h*fn1)
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)

	if ordem == 5: 
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			yn1 = y[x]+(h*((1901/720)*f[x]-(1387/360)*f[x-1]+(109/30)*f[x-2]-(637/360)*f[x-3]+(251/720)*f[x-4]))
			fn1 = solver(func, yn1, t[i])#prever com adam bashfort
			y.append((300/137)*y[x]-(300/137)*y[x-1]+(200/137)*y[x-2]-(75/137)*y[x-3]+(12/137)*y[x-4]+(60/137)*h*fn1)
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)

	if ordem == 6:
		for i in range(ordem, n+1, 1):
			x = len(y)-1
			yn1 = y[x]+(4277/1440)*h*f[x]-(2641/480)*h*f[x-1]+(4991/720)*h*f[x-2]-(3649/720)*h*f[x-3]+(959/480)*h*f[x-4]-(95/288)*h*f[x-5]
			fn1 = solver(func, yn1, t[i])#prever com adam bashfort
			y.append((360/147)*y[x]-(450/147)*y[x-1]+(400/147)*y[x-2]-(225/147)*y[x-3]+(72/147)*y[x-4]-(10/147)*y[x-5]+(60/147)*h*fn1)
			f.append(solver(func, y[x+1], t[i]))
			print(i, ' ', y[i], file=arquivo)
					
	return

def main():
	while True:
		try:
			lista = input().split(' ')
			nome_metodo = lista[0]

			if nome_metodo == 'euler':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				euler(y0, t0, h, n, func) 

			elif nome_metodo == 'euler_inverso':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				euler_inverso_prev(y0, t0, h, n, func) 		

			elif nome_metodo == 'euler_aprimorado':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				euler_aprimorado(y0, t0, h, n, func) 

			elif nome_metodo == 'runge_kutta':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				runge_kutta(y0, t0, h, n, func) 

			elif nome_metodo == 'adam_bashforth':
				x = len(lista)
				ordem = int(lista[x-1])
				y = list(range(ordem))
				for i in range(0, ordem, 1):
					y[i] = float(lista[i+1])
				t0 = float(lista[ordem+1])
				h = float(lista[ordem+2])
				n = int(lista[ordem+3])
				func = sympify(lista[ordem+4])
				adam_bashfort(y, t0, h, n, func, ordem)

			elif nome_metodo == 'adam_bashforth_by_euler':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = euler_list(y0, t0, h, ordem, func) 		
				adam_bashfort(y, t0, h, n, func, ordem)

			elif nome_metodo == 'adam_bashforth_by_euler_aprimorado':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = euler_aprimorado_list(y0, t0, h, ordem, func) 		
				adam_bashfort(y, t0, h, n, func, ordem)

			elif nome_metodo == 'adam_bashforth_by_runge_kutta':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = runge_kutta_list(y0, t0, h, ordem, func) 		
				adam_bashfort(y, t0, h, n, func, ordem)

			elif nome_metodo == 'adam_bashforth_by_euler_inverso':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = euler_inverso_list(y0, t0, h, ordem, func) 		
				adam_bashfort(y, t0, h, n, func, ordem)

			elif nome_metodo == 'adam_multon':
				x = len(lista)
				ordem = int(lista[x-1])
				y = list(range(ordem))
				for i in range(0, ordem, 1):
					y[i] = float(lista[i+1])
				t0 = float(lista[ordem+1])
				h = float(lista[ordem+2])
				n = int(lista[ordem+3])
				func = sympify(lista[ordem+4])
				adam_moulton(y, t0, h, n, func, ordem)

			elif nome_metodo == 'adam_multon_by_euler':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = euler_list(y0, t0, h, ordem, func) 		
				adam_moulton(y, t0, h, n, func, ordem)

			elif nome_metodo == 'adam_multon_by_euler_inverso':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = euler_inverso_list(y0, t0, h, ordem, func) 		
				adam_moulton(y, t0, h, n, func, ordem)

			elif nome_metodo == 'adam_multon_by_runge_kutta':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = runge_kutta_list(y0, t0, h, ordem, func) 		
				adam_moulton(y, t0, h, n, func, ordem)

			elif nome_metodo == 'adam_multon_by_euler_aprimorado':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = euler_aprimorado_list(y0, t0, h, ordem, func) 		
				adam_moulton(y, t0, h, n, func, ordem)

			elif nome_metodo == 'formula_inversa':
				x = len(lista)
				ordem = int(lista[x-1])
				y = list(range(ordem))
				for i in range(0, ordem, 1):
					y[i] = float(lista[i+1])
				t0 = float(lista[ordem+1])
				h = float(lista[ordem+2])
				n = int(lista[ordem+3])
				func = sympify(lista[ordem+4])
				formula_inversa(y, t0, h, n, func, ordem)
			
			elif nome_metodo == 'formula_inversa_by_euler':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = euler_list(y0, t0, h, ordem, func) 		
				formula_inversa(y, t0, h, n, func, ordem)

			elif nome_metodo == 'formula_inversa_by_euler_inverso':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = euler_inverso_list(y0, t0, h, ordem, func) 		
				formula_inversa(y, t0, h, n, func, ordem)

			elif nome_metodo == 'formula_inversa_by_runge_kutta':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = runge_kutta_list(y0, t0, h, ordem, func) 		
				formula_inversa(y, t0, h, n, func, ordem)

			elif nome_metodo == 'formula_inversa_by_euler_aprimorado':
				y0 = float(lista[1])
				t0 = float(lista[2])
				h = float(lista[3])
				n = int(lista[4])
				func = sympify(lista[5])
				ordem = int(lista[6])
				y = euler_aprimorado_list(y0, t0, h, ordem, func) 		
				formula_inversa(y, t0, h, n, func, ordem)

			else:
				print("try again", file=arquivo)

		except EOFError:
			break

#arquivo.close()

if __name__ == "__main__":
    main()