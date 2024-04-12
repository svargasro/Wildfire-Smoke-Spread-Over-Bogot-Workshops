# Importación de librerías necesarias
readParamsFromTable(rMean= 0.0050)
import math
from yade import pack, ymport, utils, qt
from yade.params import table
import pandas as pd



# Función para inicializar los materiales
def initialize_materials():
    """
    Define los materiales a utilizar en la simulación.
    """
    O.materials.append(FrictMat(density=rho_B, young=70e9, poisson=0.18, frictionAngle=0.1, label='spheresMat'))
    O.materials.append(FrictMat(density=2230, young=64e6, poisson=.20, frictionAngle=0.1, label='vidrio'))
    matId = O.materials.append(FrictMat())
    return matId

# Función para calcular volúmenes y masa de una esfera
def calculate_sphere_properties():
    """
    Calcula propiedades importantes de las esferas basadas en parámetros de entrada.
    """
    rMean = table.rMean  # Radio medio de las esferas
    V_o = (4/3) * math.pi * ((1000 * rMean) ** 3)  # Volumen de una esfera a escala
    V_or = (4/3) * math.pi * (rMean ** 3)  # Volumen de una esfera real
    m_o = rho_B * V_or  # Masa de una esfera
    return V_o, V_or, m_o

# Función para contar cruces y calcular flujo
def count(previousZPositions,crossings,t, flujo, k, wall,wall2, V_s, D_0):
    """
    Función ejecutada en cada iteración para contar las esferas que cruzan y verificar el balance de fuerzas.
    """
    
    flujo=0      
    print(f't={t}, flujo={flujo}, masa={flujo*m_o}, crossings={crossings}')#para poder ver ciertas cosas en la consola a la hora de modificar 
    flujo=0  
    t+=1
    k+=1
    #unbalancedforc() calcula la suma de las fuerzas desequelibradas entre las particulas
    if utils.unbalancedForce() < 0.02 and t >2:          #el t>2 es porque unbalancedforce() tiene problemas en el principio de la simulacio
        
        if wall!=None:                                   #esto se hace porque cuando intenta eliminar el muro y este no existe, hay problemas
            O.bodies.erase(wall.id)
            k = 0                         
        else:                                             # Si el muro no existe, pasar sin hacer nada
            pass

        #lo siguiente es para hallar la densidad de bulk y V_t, descomentar y cuadrar para medir
          
    if wall2 != None:
              zmax = O.bodies[wall2.id].state.pos[2]
              rint=zmax*(r/h1)
              V_t=((zmax*math.pi)/3)*(r0**2+rint**2+r0*rint)
              O.bodies[wall2.id].state.vel = Vector3(0, 0, 0)
              porosidad = yade._utils.porosity(V_t)
              Wt=35*rho_B*porosidad*(g**(1/2))*(D_0-(1.4*d))**(2.5)
              Wt2=35*rho_B*(V_s/V_t)*(g**(1/2))*(D_0-(1.4*d))**(2.5)
              print(f"zmax={zmax}, rint={rint}, V_t={V_t}, V_s={V_s}, porosidad={porosidad}, V_s/V_t={V_s/V_t}, Wt={Wt}, Wt2={Wt2}")
              
    if t==5:
        wall2 = yade.utils.wall(Vector3(0,0,90), 2, sense=-1, material=matId)
        O.bodies.append(wall2)
        O.bodies[wall2.id].state.vel = Vector3(0, 0, -1) 

     
      
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
    # Implementación de la lógica para contar cruces y actualizar el flujo y otras variables globales.

# Configuración de la simulación
def setup_simulation(matId, previousZPositions,crossings,t, flujo, k, wall,wall2, V_s, beta):
    """
    Configura los objetos, motores y condiciones iniciales de la simulación.
    """
    # Creación de esferas, muros y otros elementos de la simulación usando `pack`, `ymport`, y otros módulos de Yade.
    id_HourGl=O.bodies.append(ymport.gmsh("reloj20.mesh",scale=1000,color=(0,0,1), material='vidrio'))   
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
            [Ip2_FrictMat_FrictMat_MindlinPhys(betan=beta,betas=beta)],                          #Interacciones materiales y físicas que se considerarán en la simulación
            [Law2_ScGeom_MindlinPhys_Mindlin()],                                                 #Modelo de interacción física que se considerará en la simulación
        ),
        NewtonIntegrator(damping=0,gravity=[0,0,-g], exactAsphericalRot=True, label='newton'), #Integrador de newton que se utilizará en la simulación
        PyRunner(command='count(previousZPositions,crossings,t, flujo, k, wall,wall2, V_s)', iterPeriod=segundo, label='count', dead=False),              #Ejecuta la función count cada segundo
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


    

if __name__ == "__main__":
    
    # Inicialización de variables globales
    data = pd.DataFrame(columns=['t', 'z'])  # DataFrame para almacenar los resultados
    rho_B = 2650  # Densidad base para los materiales
    flujo = 0  # Contador de esferas que cruzan por segundo
    crossings = 0  # Contador total de esferas que cruzan
    wall = None  # Variable para manejar el muro
    k = 0  # Contador para salir del loop en equilibrio
    t = 1  # Contador de tiempo
    segundo=int(1/(0.6*(table.rMean*1000)* ((rho_B / (70 * 10**9))**(1/2))))    #cuantas iteraciones es un segundo, depende O.dt
    # segundo=int(1/O.dt)
    # print(f"segundo={segundo}")
    beta = 0.93  # Coeficiente de fricción
    wall2=None
    
    # R_0 = 1.5/100
    # R_0 = 2/100
    # R_0 = 2.5/100
    # R_0 = 3/100
    # R_0 = 3.5/100
    R_0 = 4/100
    D_0 = 2*R_0
                         
    r0 = R_0*1000
    r1 = 0.13*1000
    h1 = 0.2*1000
    r=r1-r0

    g= 9.8
    D_0 = 0.04
    r_esfera = 0.005
    d_esfera = 2*r_esfera
    
    previousZPositions = {} 
        
    V_o, V_or, m_o = calculate_sphere_properties()  # Cálculo de propiedades de las esferas
    V_s=500*V_o  
    
    matId = initialize_materials()  # Inicialización de materiales
    

    setup_simulation(matId, previousZPositions,crossings,t, flujo, k, wall,wall2, V_s, beta)  # Configuración de la simulación

    # Ejecución de la simulación
    O.run()
    if yade.utils.runningInBatch():
        # Guardado de resultados si se ejecuta en modo batch
        data.to_csv(f'Data/reloj20radio-{table.rMean}.csv', index=False)
