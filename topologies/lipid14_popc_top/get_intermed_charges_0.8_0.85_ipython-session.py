import pmx
import pmx.topology
import pmx.forcefield
ff75 = pmx.forcefield.ITPFile(fname="POPC_headgr_scaled.itp")
ff100 = pmx.forcefield.ITPFile(fname="POPC.itp")
for i,a in enumerate(ff100.atoms):
    q75 = ff75.atoms[i].q
    d = a.q - q75
    if a.q != 0.0 and d != 0.0:
        print "{id}  {atname}  {q100}  {q85}  {q80}  {q75}  {scf75}".format(id=a.id, atname=a.name, q100=a.q, q85=a.q-(d*0.6), q80=a.q-(d*0.8), q75=q75, scf75=q75/a.q)
    else:
        print "{id}  {atname}  {q100}".format(id=a.id, atname=a.name, q100=a.q)
