import odbAccess as oa
from abaqusConstants import *
import numpy as np

odb = oa.openOdb('./BRBWOPRES_IN.odb')
CT = odb.rootAssembly.instances['PART-1-1'].nodeSets['CT']
CB = odb.rootAssembly.instances['PART-1-1'].nodeSets['CB']

Toplbs = np.array([CT.nodes[i].label for i in range(len(CT.nodes))])
TopNds = np.vstack((CT.nodes[i].coordinates for i in range(len(CT.nodes))))
Botlbs = np.array([CB.nodes[i].label for i in range(len(CB.nodes))])
BotNds = np.vstack((CB.nodes[i].coordinates for i in range(len(CB.nodes))))

FEs = {"FACE1": list([0, 1, 2, 3]),
       "FACE2": list([4, 7, 6, 5]),
       "FACE3": list([0, 4, 5, 1]),
       "FACE4": list([1, 5, 6, 2]),
       "FACE5": list([2, 6, 7, 3]),
       "FACE6": list([3, 7, 4, 0])}
INTELEMCVT = []
for i in range(1, 6):
    CTs = odb.rootAssembly.surfaces['CT%dS' % (i)]
    [INTELEMCVT.extend([list(np.array(CTs.elements[0][i].connectivity)[FEs[CTs.faces[0][i].getText()]])]) for i in range(len(CTs.elements[0]))]
INTELEMCVT = np.array(INTELEMCVT)
Ne = INTELEMCVT.shape[0]
INTELS = np.zeros((Ne, 5), dtype="int")
for i in range(Ne):  # Extracting Element Node Indices in C Indexing (0->N-1)
    INTELS[i, 0] = i
    INTELS[i, 1:] = [np.where(Toplbs == INTELEMCVT[i][j])[0][0] for j in range(3, -1, -1)]
INTELS += 1  # Conversion to Fortran Indexing (1->N)

# Text file Output
np.savetxt('Nodes.dat', TopNds[:, 0:2])
np.savetxt('Elements.dat', INTELS, fmt = "%d")
