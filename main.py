#METODOS NUMERICOS - Act. 6
#Sistema de ecuaciones de la forma:
# ax + by + cz = d
from calculator import Gui
#from activity2 import guiTayler
#from activity3 import guiRoots
#from activity4 import guiPolyInter
from activity6 import guiSysLinEq

api = Gui('Actividad 6', 'METODOS NUMERICOS', 'gray')
#api.addActivity('Series de Taylor', guiTayler, '#009ece')
#api.addActivity('Calculo de Raices', guiRoots, '#009ece')
#api.addActivity('Interpolacion Polinomica', guiPolyInter, '#009ece')
api.addActivity('Iniciar', guiSysLinEq, '#009ece')
api.mainloop()