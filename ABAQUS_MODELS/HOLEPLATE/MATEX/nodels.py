from abaqus import *
from abaqusConstants import *
import regionToolset
import job
import step
import sets

import numpy as np

mdb = openMdb('../HOLEPLATE.cae')
mdl = mdb.models['Model-1']
ras = mdl.rootAssembly

Ss = mdl.parts['Part-1'].surfaces
Vs = Ss.keys()

# Node, Element IDs
Ndi = []
Eli = []
Elti = []
for s in Vs:
    for i in range(len(Ss[s].nodes)):
        Ndi.append(Ss[s].nodes[i].label)
    Eli.append([Ss[s].elements[i].label for i in range(len(Ss[s].elements))])
    for ne in Eli[-1]:
        Elti.append(ne)
Ndi = np.unique(Ndi)
Elti = np.unique(Elti)

# Nodal Locations
Nds = np.array([mdl.parts['Part-1'].nodes[n-1].coordinates for n in Ndi])

# Node Dictionary
Nd_dict = dict(zip(Ndi, range(len(Ndi))))
El_dict = dict(zip(Elti, range(len(Elti))))

# Element Connectivities
Elcs = np.zeros((len(Elti), 6), 'int')
ke = 0
for ne in Elti:
    fi = np.argwhere([all([(x in Ndi) for x in [mdl.parts['Part-1'].elements[ne-1].getElemFaces()[f].getNodes()[i].label for i in range(4)]]) for f in range(6)])[0][0]

    Elcs[ke, 0] = ke
    Elcs[ke, 1:-1] = [Nd_dict[x] for x in [mdl.parts['Part-1'].elements[ne-1].getElemFaces()[fi].getNodes()[i].label for i in range(4)]]
    Elcs[ke, -1] = np.argwhere([ne in Eli[i] for i in range(5)])[0][0]
    ke += 1

    print('Done %d' % (ke))

# Print to Text
np.savetxt('Nodes.dat', Nds)
np.savetxt('Elements.dat', Elcs)
np.savetxt('NodeLabels.dat', Ndi)
