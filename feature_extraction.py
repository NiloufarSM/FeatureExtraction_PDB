import sys
from Bio.PDB import *
import csv

DEBUG = False
if DEBUG :
    pdbfile = "1hl3.pdb"
else:
    arguments = sys.argv
    pdbfile = arguments[1]
    #print(arguments)
pdbname= pdbfile.replace('.pdb','')
with open(pdbname+'.csv', 'w', newline='') as out:
    output = csv.writer(out, delimiter=' ')
    #out.write(pdbfile)
    if DEBUG: print(pdbname)
    out.write("Residue\tChain\tSS_DSSP\tEXP_HSE_B_U\tEXP_HSE_B_D\tEXP_DSSP_RASA")
    output.writerow("")

    p = PDBParser()
    structure = p.get_structure("one", pdbfile)
    model = structure[0]
    for chain in model:
        #chain = model["A"]
        if DEBUG: print(len(chain))
        dssp = DSSP(model, pdbfile,dssp='mkdssp')

        # forth item is relative accessible surface area
        #for item in list(dssp.keys()):
        #    print(dssp[item][3])


        #calculate HSE based on the approximate CA-CB vectors.
        #hsea = HSExposureCA(model,radius=12)
        #calculate HSE based on the real CA-CB vectors.
        hseb = HSExposureCB(model,radius=12)
        #Residue exposure as number of CA atoms around its CA atom.
        #hsen = HSExposure.ExposureCN(model,radius=12)

        residues = chain.get_residues()
        if DEBUG:print(residues)
        for c, r in enumerate(residues):
            # The first and the last residues do not have HSE measure
            """if c == 1:
                if DEBUG: print("skip first one")
                c = c + 1
                continue
            if c == len(chain):
                if DEBUG: print("skip last one")
                break"""
            c = c + 1
            try:
                if is_aa(r):
                    if DEBUG: print("residue", r.get_resname())
                    if DEBUG: print(r.xtra)
                    #if DEBUG: print("EXP_HSE_B_U:", r.xtra["EXP_HSE_B_U"] , "EXP_HSE_B_D:" , r.xtra["EXP_HSE_B_D"] , "EXP_DSSP_RASA:",r.xtra["EXP_DSSP_RASA"] )
                    out.write(r.get_resname()+"\t"+chain.get_id()+"\t")
                    try:
                        out.write(str(r.xtra["SS_DSSP"])+"\t")
                    except:
                        out.write("NA\t")
                    try:
                        out.write(str(r.xtra["EXP_HSE_B_U"])+"\t")
                    except:
                        out.write("NA\t")
                    try:
                        out.write(str(r.xtra["EXP_HSE_B_D"])+"\t")
                    except:
                        out.write("NA\t")
                    try:
                        out.write(str(r.xtra["EXP_DSSP_RASA"]))
                    except:
                        out.write("NA")

                    output.writerow("")
                    #print("A:", hsea[(chain.get_id(), r.get_id())])
                    #print("B:", hseb[(chain.get_id(), r.get_id())])
                    #print("N:", hsen[(chain.get_id(), r.get_id())])

                    # print(hsea[(chain.get_id(),r.get_id())][2])
                    # print(hseb[(chain.get_id(),r.get_id())][2])

                    """
                    try:
                        #contact number
                        print(r.xtra["EXP_CN"])
                        #HSE alpha up
                        print(r.xtra["EXP_HSE_A_U"])
                        print()
                    except:
                        pass
                    """
            except Exception as ex:
                template = "An exception of type {0} occurred. Arguments:\n{1!r}"
                message = template.format(type(ex).__name__, ex.args)
                print (message)
                pass
