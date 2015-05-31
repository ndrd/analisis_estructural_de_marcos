from numpy.linalg import inv
from numpy import transpose
from numpy import concatenate as concat
import numpy as np
import sys
import math
from pprint import pprint as pp

nodos = []
barras = []
K_matriz = None
R_matriz = []

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
		self.fx = _datos[4]
		self.fy = _datos[5]
		self.momento = _datos[6]

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

		self.crear_matriz_transformacion()
		self.crear_matriz_desconocida()
		self.crear_matriz_K_rigidez()

	def crear_matriz_transformacion():

		puntoA = self.puntoA
		puntoB = self.puntoB

		L = math.sqrt( math.pow(puntoA.x - puntoB.x, 2) + math.pow(puntoA.y - puntoB.y, 2) )
		S = ( puntoB.y - puntoA.y ) / L
		C = ( puntoB.x - puntoA.x ) / L

		self.matriz_transformada = [
			[ C,-S, 0],
			[ S, C, 0],
			[ 0, 0, 1]
		]

		self.matriz_transformada_transpuesta =  transpose(self.matriz_transformada)
		self.L = L

	def crear_matriz_desconocida():

		
