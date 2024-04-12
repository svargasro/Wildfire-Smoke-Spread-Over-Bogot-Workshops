# yade-batch -j 2 --job-threads=1 parametros.txt ultimaprueba.py
readParamsFromTable(rMean= 0.0050)
from yade.params import table 

import pandas as pd
data = pd.DataFrame(columns=['t', 'z'])

rho_B = 2650

#definir materilaes y su propiedades
O.materials.append(FrictMat(density=rho_B,young=70e9,poisson=0.18,frictionAngle=0.1, label='spheresMat'))
O.materials.append(FrictMat(density=2230,young=64e6, poisson=.20,frictionAngle=0.1,label='vidrio'))
#O.materials.append(FrictMat(density=0.00265,young=(70e9),poisson=0.18,frictionAngle=0.1, label='spheresMat'))
#O.materials.append(FrictMat(density=0.002230,young=(64e6), poisson=.20,frictionAngle=0.1,label='vidrio'))
matId = O.materials.append(FrictMat())

#variables que se utilizaran, y algunas que se necesitan inicializar antes de ser llamadas

V_o=(4/3)*math.pi*((1000*table.rMean)**3)         #volumen de una esfera a escala
V_or=(4/3)*math.pi*((table.rMean)**3)             #volumen de una esfera real
m_o=rho_B*V_or                                      #masa de una esfera

flujo=0                  #cantidad de esferas que pasan en un segundo de un lado al otro del reloj   
crossings=0              #cantidad total de esferas que cruzan no es mas que la suma del flujo, se hizo para confirmar que al final de la simulacion todas las esferas pasaran

wall = None              #se inicializa el muro 
k=0                      #variable que sirve para salir del programa la segunda vez que las esferas queden casi en equilibrio
t=1                      #contador de tiempo
segundo=int(1/(0.6*(table.rMean*1000)* ((rho_B / (70 * 10**9))**(1/2))))    #cuantas iteraciones es un segundo, depende O.dt
# segundo=int(1/O.dt)
# print(f"segundo={segundo}")


#las siguientes variables se activan cuando se quiera medir la densidad de bulk

wall2=None
V_s=500*V_o
# R_0 = 1.5/100
# R_0 = 2/100
# R_0 = 2.5/100
# R_0 = 3/100
# R_0 = 3.5/100
R_0 = 4/100
D_0 = 2*R_0

# r = 0.55/100
# d = 2*r

# d = 0.011
# d = 0.006
# d = 0.008
# D_0=0.04
d = 0.011                                       
r0 = R_0*1000
r1 = 0.13*1000
h1 = 0.2*1000
r=r1-r0

g= 9.8


print(f"D_0={D_0}, d={d}")
#para almacenar las posiciones en z, en el intervalo anterior

previousZPositions = {}


#funcion para contar cuantas esferas han cruzado y verificar el balance de fuerzas
def count():  
    global previousZPositions
    global crossings,t, flujo, V_or, k, wall, m_o,wall2, V_s, R_0, r0, r1, h1, r, g, D_0, d 
    t+=1
    k+=1
    t_estable=2 #Tiempo en el que se considera que las esferas estan en equilibrio
    print(f't={t}, flujo={flujo}, masa={flujo*m_o}, crossings={crossings}')#para poder ver ciertas cosas en la consola a la hora de modificar
    #unbalancedforc() calcula la suma de las fuerzas desequelibradas entre las particulas
    if utils.unbalancedForce() < 0.02 and t >t_estable:          #el t>2 es porque unbalancedforce() tiene problemas en el principio de la simulacio
        
        if wall!=None:                                   #esto se hace porque cuando intenta eliminar el muro y este no existe, hay problemas
            O.bodies.erase(wall.id)
            k = 0                         
        else:                                             # Si el muro no existe, pasar sin hacer nada
            pass

        #lo siguiente es para hallar la densidad de bulk y V_t, descomentar y cuadrar para medir
          
    if  t==3:
        wall2 = yade.utils.wall(Vector3(0,0,90), 2, sense=-1, material=matId)
        O.bodies.append(wall2)
        O.bodies[wall2.id].state.vel = Vector3(0, 0, -1) 
        densidadbulk()

        

     
      
    #cuando se vuelve a quedar en equilibrio en el fondo del reloj
    if k > 5000 and utils.unbalancedForce() <  0.02:                  
        O.pause()
    currentZPositions = {b.id: b.state.pos[2] for b in O.bodies if isinstance(b.shape, Sphere)}
    #print(crossings)
    # Itera sobre las esferas y verifica los cambios de positivo a negativo en la posición z
    for bodyId, currentZ in currentZPositions.items():
        if bodyId in previousZPositions:
            previousZ = previousZPositions[bodyId]
            if previousZ > 0 and currentZ <= 0:
                flujo += 1
    
    crossings+=flujo

    # Actualiza las posiciones z anteriores
    previousZPositions = currentZPositions
 
    global data
    if yade.utils.runningInBatch():
        #z = O.bodies[sphereId].state.pos[2]
        data.loc[len(data)] = [O.dt*O.iter, flujo*m_o]

def densidadbulk():
        zmax = O.bodies[wall2.id].state.pos[2]
        rint=zmax*(r/h1)
        V_t=((zmax*math.pi)/3)*(r0**2+rint**2+r0*rint)
        O.bodies[wall2.id].state.vel = Vector3(0, 0, 0)
        porosidad = yade._utils.porosity(V_t)
        Wt=35*rho_B*porosidad*(g**(1/2))*(D_0-(1.4*d))**(2.5)
        Wt2=35*rho_B*(V_s/V_t)*(g**(1/2))*(D_0-(1.4*d))**(2.5)
    #   print(f"zmax={zmax}, rint={rint}, V_t={V_t}, V_s={V_s}, porosidad={porosidad}, V_s/V_t={V_s/V_t}, Wt={Wt}, Wt2={Wt2}")
        print(f"zmax={zmax}, r0={r0}, rint={rint}, porosidad={porosidad}")

from yade import ymport
#crear el reloj
id_HourGl=O.bodies.append(ymport.gmsh("reloj" + str(int(r0)) + ".mesh",scale=1000,color=(0,0,1), material='vidrio'))   
#generar un paquete de esferas
from yade import pack
sp= pack.SpherePack()
sp.makeCloud((-65,-65,144),(65,65,200), num=400, rMean=table.rMean*1000, rRelFuzz=0.0, seed=2)
sp.toSimulation(material='spheresMat')
#generar otro paquete de esferas
sp1= pack.SpherePack()
sp1.makeCloud((-50,-50,80),(50,50,140), num=100,rMean=table.rMean*1000, rRelFuzz=0.0, seed=2)
sp1.toSimulation(material='spheresMat')


#crear el muro
wall = yade.utils.wall(Vector3(0,0,0), 2, sense=0, material = matId)
O.bodies.append(wall)




#motores y demas elementos de la simulaciom

O.engines=[
    ForceResetter(),                                                                         #Reinicia las fuerzas en cada paso de tiempo
    InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb(),Bo1_Wall_Aabb()]),             #Determina las colisiones entre los diferentes tipos de objetos
    InteractionLoop(                                                                         #Establece las interacciones entre diferentes tipos de objetos en la simulación
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(), Ig2_Wall_Sphere_ScGeom()],   #Interacciones geométricas que se considerarán en la simulación
        [Ip2_FrictMat_FrictMat_MindlinPhys(betan=0.93,betas=0.93)],                          #Interacciones materiales y físicas que se considerarán en la simulación
        [Law2_ScGeom_MindlinPhys_Mindlin()],                                                 #Modelo de interacción física que se considerará en la simulación
    ),
    NewtonIntegrator(damping=0,gravity=[0,0,-g], exactAsphericalRot=True, label='newton'), #Integrador de newton que se utilizará en la simulación
    PyRunner(command='count()', iterPeriod=segundo, label='count', dead=False),              #Ejecuta la función count cada segundo
    #PyRunner(command='densidadbulk()', iterPeriod=10, label='densidad', dead=False),

]
#importar el modulo qt
#inicializacion de los objeros renderer y View de qt por si no se corre yade-batch
if not yade.utils.runningInBatch():
    try:
        from yade import qt
        rndr = yade.qt.Renderer()
        v = qt.View()
    except: pass

O.dt=0.6*PWaveTimeStep()
O.stopAtTime = 105000 * O.dt
O.run()
waitIfBatch()

if yade.utils.runningInBatch():
    data.to_csv(f'Data/reloj20radio-{table.rMean}.csv', index=False)
