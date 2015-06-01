from numpy.linalg import inv
from numpy import transpose
from numpy import concatenate as concat
import numpy as np
import sys
import math
from pprint import pprint as pp
import json

nodos = []
barras = {}
K_matriz = None
R_matriz = []

identidad_negativa = [
	[-1,0,0],
	[0,-1,0],
	[0,0,-1]
]

def buscar_nodo_por_id(identificador):
	for _nodo in nodos:
		if _nodo.id == identificador:
			return _nodo

class Nodo(object):
	"""Abstraccion para trabajar con nodos"""
	
	def __init__(self, _datos):
		self.id = _datos[1]
		self.x  = float(_datos[2])
		self.y  = float(_datos[3])

	def extraer_fuerzas( self, R ):

		inicio = (int(self.id) - 1 ) * 3

		#Si la fuerza esta definida
		if self.fx != 'u':
			R[inicio] = float(self.fx)
		else:
			R[inicio] = 'R' + str(inicio)

		#Si la fuerza esta definida
		if self.fy != 'u':
			R[inicio+1] = float(self.fy)
		else:
			R[inicio+1] = 'R' + str(inicio)

		if self.momento != 'u':
			R[inicio+2] = float(self.momento)
		else:
			R[inicio+2] = 'R'+str(inicio+2)


	def extraer_fuerzas_r( self, r ):

		inicio = (int(self.id) - 1 ) * 3

		#Si la fuerza esta definida
		if self.fx != 'u':
			r[inicio] = float(self.fx)
		else:
			r[inicio] = 'R' + str(inicio)

		#Si la fuerza esta definida
		if self.fy != 'u':
			r[inicio+1] = float(self.fy)
		else:
			r[inicio+1] = 'R' + str(inicio)

		if self.momento != 'u':
		  	r[inicio+2] = float(self.momento)
		else:
		  	r[inicio+2] = 'R'+str(inicio+2)

class Barra(object):

	""" abstraccion de la barra """

	def __init__(self, _datos):

		self.puntoA = buscar_nodo_por_id(_datos[1])
		self.tipoPuntoA = _datos[2]

		self.puntoB = buscar_nodo_por_id(_datos[3])
		self.tipoPuntoB = _datos[4]

		self.fEmpotramientoA  =  0	
		self.fEmpotramientoB  =  0	

		self.FxA = 0
		self.FxB = 0

		self.FyA = 0
		self.FyB = 0
 
		self.E = float(_datos[5])
		self.I = float(_datos[6])
		self.A = float(_datos[7])

		self.cargas_distribuidas = float(_datos[8])

		self.FxLocal = float(_datos[9])
		self.FyLocal = float(_datos[10])
		self.distancia_fuerza_externa = float(_datos[11])

		self.id  = int(_datos[12])

		self.crear_matriz_transformacion()
		self.crear_matriz_rigidez_elemento()
		self.crear_matriz_K_ceros()


		# las cargas distribuidas son cero 
		if self.cargas_distribuidas != "u":
			self.calcular_fuerzas_distribuidas()
		else:
			self.calcular_fuerzas_empotramiento_en_y()




	def crear_matriz_transformacion(self):

		puntoA = self.puntoA
		puntoB = self.puntoB

		L = math.sqrt( math.pow(puntoA.x - puntoB.x, 2) + math.pow(puntoA.y - puntoB.y, 2) )
		S = ( puntoB.y - puntoA.y ) / L
		C = ( puntoB.x - puntoA.x ) / L

		self.matriz_transformada = [
			[ C, S, 0],
			[-S, C, 0],
			[ 0, 0, 1]
		]

		self.matriz_transformada = np.matrix(self.matriz_transformada)
		self.matriz_transformada_transpuesta =  self.matriz_transformada.transpose()

		self.L = L

	def crear_matriz_rigidez_elemento(self):
		#Matriz de fuerza de cada nudo
		self.Km = []
		EA = self.E * self.A
		EI = self.E * self.I
		L = self.L

		self.KAA = np.matrix([
				[  EA/L    ,      0      ,     0     ],
				[    0     , 12*EI/L**3  , 6*EI/L**2  ], 
				[    0     , 6*EI/L**2   , 4*EI/L      ]
		])
		self.KAAI = np.dot(np.dot( self.matriz_transformada_transpuesta, self.KAA), self.matriz_transformada)

		self.KAB = np.matrix([
				[ -EA/L    ,      0      ,     0      ],
				[   0      , -12*EI/L**3 ,  6*EI/L**2 ],
				[   0      , -6*EI/L**2  ,  2*EI/L    ],
		])
		self.KABI = np.dot(np.dot( self.matriz_transformada_transpuesta, self.KAB), self.matriz_transformada)

		self.KBA = np.matrix([
				[  -EA/L   ,      0      ,     0     ]  ,
				[    0     , -12*EI/L**3 , -6*EI/L**2 ] ,
				[    0     , 6*EI/L**2   , 2*EI/L     ]
		])
		self.KBAI = np.dot(np.dot( self.matriz_transformada_transpuesta, self.KBA), self.matriz_transformada)

		self.KBB = np.matrix([
				[ EA/L    ,      0      ,     0     ],
				[ 0      , 12*EI/L**3  , -6*EI/L**2 ],
				[ 0      , -6*EI/L**2  ,  4*EI/L    ]
		])
		self.KBBI = np.dot(np.dot( self.matriz_transformada_transpuesta, self.KBB), self.matriz_transformada)

		



	def crear_matriz_K_ceros(self):

		ceros = [
			[0,0,0],
			[0,0,0],
			[0,0,0],
		]

	def calcular_fuerzas_empotramiento_en_y(self):
		
		a = self.distancia_fuerza_externa
		b = self.L - a
		P = self.FyLocal 
		L = self.L

		Fx = -self.FxLocal / 2

		M1 = ( P *  a * b ** 2 ) / L ** 2
		M2 = (-P *  a ** 2 * b ) / L ** 2

		Fy1 = ((P * b ** 2 ) / L ** 3 ) * (3 * a + b)
		Fy2 = ((P * a ** 2 ) / L ** 3 ) * (3 * b + a)

		self.matriz_fuerzas_empotramiento_A = [[Fx], [Fy1], [M1]]
		self.matriz_fuerzas_empotramiento_B = [[Fx], [Fy2], [M2]]

	def calcular_fuerzas_distribuidas(self):
				
		W = self.cargas_distribuidas 
		L = self.L

		M1 = ( W * L ** 2 ) / 12
		M2 = -M1

		Fy1 = (W * L) /  2
		Fy2 = Fy1

		self.matriz_fuerzas_empotramiento_A = [[0], [Fy1], [M1]]
		self.matriz_fuerzas_empotramiento_B = [[0], [Fy2], [M2]]



def mostrar_nodos():
	for nodo in nodos:
		print json.dumps(nodo.__dict__)

def mostrar_barras():
	for barra in barras:
		print json.dumps(barra.__dict__)




def main(argv):
	#leemos los datos del archivo de entrada
	datos = [linea.rstrip('\n').split() for linea in open(argv[1])]
	#para cada renglon construimos la barra o el nodo, segun sea el caso
	for _dato in datos:

		if len(_dato) != 0:
			if _dato[0] == 'nodo':
				nodos.append(Nodo(_dato))
			elif _dato[0] == 'barra':
				barra = Barra(_dato)
				barras.update({barra.id : barra})
	#se determina el numero de barras que comparten un nodo			
	barras_por_nodo = {}
		
	for barra in barras.values():
		if barra.puntoA.id in barras_por_nodo:
			barras_por_nodo[barra.puntoA.id] += 1
		else:
			barras_por_nodo.update({barra.puntoA.id : 1})
		if barra.puntoB.id in barras_por_nodo:
			barras_por_nodo[barra.puntoB.id] += 1
		else:
			barras_por_nodo.update({barra.puntoB.id : 1})

	#se determina cuales son los nodos que comparten 2 o mas barras	
	nodos_comunes = []

	for nodo in barras_por_nodo.keys():
		if barras_por_nodo[nodo] >= 2:
			nodos_comunes.append(nodo)

	#se busca el extremo de la barra que participa		
	nodos_a_procesar = {}
	nodos_comunes.sort()

	for nodo in nodos_comunes:

		lado_barra = {}
		#se selecciona si es el lado A o B de la barra (el que participa)	
		for barra in barras.values():
			if barra.puntoA.id == nodo:
				lado_barra.update({barra.id : 'a' })
			if barra.puntoB.id == nodo:
				lado_barra.update({barra.id : 'b' })


	#los nodos a procesar seran los que usaremos para encontrar las Ffij			
		if nodo in nodos_a_procesar.keys():
			nodos_a_procesar[nodo].update({lado_barra})
		else:
			nodos_a_procesar[nodo] = lado_barra


	#se determinara la matriz de fuerzas de fijacion para cada nodo de interes 
	#(posteriormente se uniran las ncesarias para formar el vector Ffij)	
	matrices_de_fijacion = []
	matriz_final_de_fijacion = None
	matriz_fuerzas_efectivas = None


	print json.dumps(nodos_a_procesar, indent=4, sort_keys=True)

	nodos_a =  nodos_a_procesar.keys()
	nodos_a.sort()


	for nodo in nodos_a:

		lado_barra = nodos_a_procesar[nodo]

		matrices_suma  = np.zeros((3,1))

		for id_barra in lado_barra.keys():
			#determinamos el lado de la barra a operar
			if lado_barra[id_barra] == 'a':
				barra = barras[id_barra]
				matrizFE = barra.matriz_fuerzas_empotramiento_A
			else:
				barra = barras[id_barra]
				matrizFE = barra.matriz_fuerzas_empotramiento_B

			matriz_transformada =  np.matrix(barra.matriz_transformada)
			matriz_transformada_t =  matriz_transformada.transpose()
			matrizFE = (np.matrix(matrizFE))

			m_temp = np.dot(matriz_transformada_t, matrizFE)
			matrices_suma =  matrices_suma + m_temp

		matrices_de_fijacion.append(matrices_suma)

	# unimos las matrices de fijacion en una coluna
	matriz_final_de_fijacion = np.concatenate((tuple(matrices_de_fijacion)), axis=0)
	matriz_fuerzas_efectivas = matriz_final_de_fijacion * -1

	for i in barras.keys():
		barra = barras[i]
		barra.crear_matriz_rigidez_elemento()

if __name__ == '__main__':
	main(sys.argv)		
