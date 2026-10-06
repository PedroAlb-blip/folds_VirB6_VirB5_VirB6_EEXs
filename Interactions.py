import matplotlib.pyplot as plt
import os
import re

resdict = {"ALA": "A","CYS": "C","ASP": "D","GLU": "E","PHE": "F","GLY": "G","HIS": "H","ILE": "I","LYS": "K","LEU": "L","MET": "M","ASN": "N","PRO": "P","GLN": "Q","ARG": "R","SER": "S","THR": "T","VAL": "V","TRP": "W","TYR": "Y"}
pMBUI4EEX = "CKTEETMTVDYYLANPDARKQKLVQCANNPGDEAIQANCVNAKEAQHKSLFDPDNTGMPSIK"
R388EEX = "CKDEPTHDVQYYLDNPKERTEKIAECKNNPGEKALAANCVNATEAQRKAMFQGDGMPKIR"

pMBUIVirB5_aligned = "-QIPVTDGASIAQNMANQVQTMAQWAQQLMEMQQQYTELQRHWEKMDEMEKLVTGARYLGVP-NTSGVNGQIPGTVQD--------IYSGSYGTTQGIMDAERINGTTATTQDALNERYLKNAAAEKAAAINTYQGAQERLANINTLISKIDSAQDPKAIQDLQARIAAEQALVQNETTKLQMVSQLASAEQKLVEDQRRQMNRDILNPDNTGMPEIK"
pMBUIVirB6 = "MAFSLFAPLFEKVDTATQAFVSDISANAIATITPFVTLGLTLGFIAYAILIIRGAVDMPILDFLARSLRIGIVVSVALAGGLYQSQIASAIIELPNSLAVALITDPTEGAGAASIIDVAAGKGFDAAGRAFQQSSFFSTDGLLYGIFGIMIILATAVIVAIGGAFILLAKLALALLAGLGPLFIVGLLWQPTSRFFEMWVAQVLNYVLLVVLFSAVFGLMMSIYGGYVEGVKFESGVNVAYSLGGAFILSAAMVLILMQLPGIAGALAGGVGISYMWELRSLRGGASSVGRGLAGAGRSIYGRGAKGADGSRGAATGMAPAAMRGAQKAYGYFKGSKAA"
R388VirB5_aligned = "QGIPVFDGTRALD-----------FVQQFARMKEQLDTAKDQLAEAQRMYEAVTGGRGLGDLMRNAQLREYLPDDLRTVYDSANGGGYSGISGSINDILRDERLNGSVADMRRSIEERSRTARTTDKAVGCGPTKGHSN---ASRRLKDEISRTQDQKAIEELQARIAGEQAAIQNETTKLQMIAQLRQAEQALISEQRRERNMPILSSGNQGMPTIQ"
R388VirB6 = "MAFELFTPLFNKIDQTTATYVTDISSRAIAAITPVVSVGLTLGFITYGWLIIRGAVEMPVAEFLNRCLRIGIIVSIALAGGLYQGEIANAITTVPDELASALLGNPTEGASAAALVDQSAQQGFDRASEAFEEAGFFSSDGLLYGLFGIIIPAGYGPAGAIGGAFLLLAKIALALLAGLGPLFILALIWQPTHRFFDQWAQQVLNYGLLIVLFAAVFGLLMQIFGSYMADLRFDGAQNVAYAIGGSVILSIVSIVLLMQLPSIASGLAGGIGLGYMWELRSMRSGAGAAMRGGRAMARGARAAPGAARGAAVGAANMAKTVATGGAGVARAAAGYFRGRKAG"

pMBUIVirB5 = "QIPVTDGASIAQNMANQVQTMAQWAQQLMEMQQQYTELQRHWEKMDEMEKLVTGARYLGVPNTSGVNGQIPGTVQDIYSGSYGTTQGIMDAERINGTTATTQDALNERYLKNAAAEKAAAINTYQGAQERLANINTLISKIDSAQDPKAIQDLQARIAAEQALVQNETTKLQMVSQLASAEQKLVEDQRRQMNRDILNPDNTGMPEIK"
pMBUIVirB6_aligned = "MAFSLFAPLFEKVDTATQAFVSDISANAIATITPFVTLGLTLGFIAYAILIIRGAVDMPILDFLARSLRIGIVVSVALAGGLYQSQIASAIIELPNSLAVALITDPTEGAGAASIIDVAAGKGFDAAGRAFQQSSFFSTDGLLYGIFGIMIILATAVIVAIGGAFILLAKLALALLAGLGPLFIVGLLWQPTSRFFEMWVAQVLNYVLLVVLFSAVFGLMMSIYGGYVEGVKFESGVNVAYSLGGAFILSAAMVLILMQLPGIAGALAGGVGISYMWELRSLRGGASSVGRGLAGAGRSIYGRGAKGADG-SRGAATGM-------APAAMRGAQKAYGYFKGSKAA"
R388VirB5 = "QGIPVFDGTRALDFVQQFARMKEQLDTAKDQLAEAQRMYEAVTGGRGLGDLMRNAQLREYLPDDLRTVYDSANGGGYSGISGSINDILRDERLNGSVADMRRSIEERSRTARTTDKAVGCGPTKGHSNASRRLKDEISRTQDQKAIEELQARIAGEQAAIQNETTKLQMIAQLRQAEQALISEQRRERNMPILSSGNQGMPTIQ"
R388VirB6_aligned = "MAFELFTPLFNKIDQTTATYVTDISSRAIAAITPVVSVGLTLGFITYGWLIIRGAVEMPVAEFLNRCLRIGIIVSIALAGGLYQGEIANAITTVPDELASALLGNPTEGASAAALVDQSAQQGFDRASEAFEEAGFFSSDGLLYGLFGIIIPAGYGPAGAIGGAFLLLAKIALALLAGLGPLFILALIWQPTHRFFDQWAQQVLNYGLLIVLFAAVFGLLMQIFGSYMADLRFDGAQNVAYAIGGSVILSIVSIVLLMQLPSIASGLAGGIGLGYMWELRSMRSGAGAAMRGGRAM-----ARGARAAPGAARGAAVGAANMAKTVATGGAGVARAAAGYFRGRKAG"


Vir_list = [pMBUIVirB5, pMBUIVirB6, R388VirB5, R388VirB6]
Vir_list_aligned = [pMBUIVirB5_aligned, pMBUIVirB6_aligned, R388VirB5_aligned, R388VirB6_aligned]
loi = []
interactome = {"EEX": [], "Vir": []}

for i in os.listdir():
    if os.path.isdir(i):
        if "Interactions_PDBePISA.txt" in os.listdir(i):
            file = open(os.path.join(i, "Interactions_PDBePISA.txt"), 'r')
            loi.append(i)
            interactions = file.readlines()
            counter = 0
            othercounter = 0
            for j in range(0,len(interactions)):
                if "CHAIN" in interactions[j]:
                    EEX = interactions[j].split(" ")[1]
                    Vir = interactions[j].split(" ")[2].split("\n")[0]
                    interactome["EEX"].append(EEX)
                    interactome["Vir"].append(Vir)
                    loi.append(EEX)
                elif "SALTBRIDGES" in interactions[j]:
                    index_2 = j
                    counter += 1
                    if counter%2 == 0:
                        loi.append("SALTBRIDGES")
                        loi.append(interactions[index_1:index_2])
                    index_1 = j
                elif "HYDROGENBONDS" in interactions[j]:
                    index_3 = j
                    othercounter += 1
                    if othercounter%2 == 0:
                        loi.append("HYDROGENBONDS")
                        loi.append(interactions[index_4:index_3])
                    index_4 = j

loi.append("END")
fullist = []
onelist = []
for i in range(0, len(loi)):
    if type(loi[i]) == list:
        onelist=onelist + loi[i]
    elif type(loi[i]) == str and loi[i].startswith("2026") and onelist != []:
        fullist.append(onelist)
        onelist = []
    elif i == len(loi) - 1:
        fullist.append(onelist)



results = open("interacting_residues.txt", 'w')
pMBUI_residues = []
R388_residues = []
iteration = 0

for i in loi:
    if type(i) == str and len(i) == 1: # in case its a single letter belonging or refering to the chain of the eex protein
        key = i
        #interacting = interactome(iteration)
        iteration += 1
        results.write(key + "\n")
    elif type(i) == str and len(i) > 1:
        results.write(i + "\n")
        if "eex_pmbui" in i:
            eex = "pMBUI"
        elif "eex_r388" in i:
            eex = "R388"
    elif type(i) == list:
        pos = 0
        residues = []
        interacting_prot = []
        for j in i:
            #print(pos)
            if key + ":" in j and "[ O ]" not in j and "[ N ]" not in j and "[" in j:
                res = re.findall(r' [0-9]+\[', j)
                if "STRUCTURE2" in j:
                    results.write(i[pos-2])
                    interacting_prot.append(i[pos-2])
                    results.write(i[pos-1])
                    results.write(i[pos])
                    results.write(res[0].replace(" ", "").split("[")[0] + "\n")
                    residues.append(int(res[0].replace(" ", "").split("[")[0]))
                elif "STRUCTURE1" in j:
                    results.write(res[0].replace(" ", "").split("[")[0] + "\n")
                    results.write(i[pos])
                    results.write(i[pos+1])
                    interacting_prot.append(i[pos+2])
                    results.write(i[pos+2])
                    residues.append(int(res[0].replace(" ", "").split("[")[0]))
            pos += 1
        prot =0
        
        for k in residues:
            if eex == "pMBUI":
                pMBUI4EEX = pMBUI4EEX[0 : k - 1] + pMBUI4EEX[k - 1].lower() + pMBUI4EEX[k : len(pMBUI4EEX)]
                prot = 1
            elif eex == "R388":
                R388EEX = R388EEX[0 : k - 1] + R388EEX[k - 1].lower() + R388EEX[k : len(R388EEX)]
                prot = 0
        if prot == 1:
            results.write(pMBUI4EEX + "\n")
            pMBUI_residues.append(pMBUI4EEX)
            pMBUI4EEX = pMBUI4EEX.upper()
        elif prot == 0:
            results.write(R388EEX + "\n")
            R388_residues.append(R388EEX)
            R388EEX = R388EEX.upper()



pMBUI_heat = [0]*len(pMBUI_residues[0])
R388_heat = [0]*len(R388_residues[0])
for i in range(0, len(pMBUI_residues)):
    index = 0
    for j in pMBUI_residues[i]:
        if j.islower():
            pMBUI_heat[index] = pMBUI_heat[index] + 1
        index = index + 1

for i in range(0, len(R388_residues)):
    index = 0
    for j in R388_residues[i]:
            if j.islower():
                R388_heat[index] = R388_heat[index] + 1
            index = index + 1


results.write("\n")
results.write("pMBUI4 HEATMAP OF INTERACTIONS")
results.write("\n")
results.write(str(pMBUI_heat))
results.write("\n")
results.write("R388 HEATMAP OF INTERACTIONS")
results.write("\n")
results.write(str(R388_heat))

R388EEX = "CKDEPTHDVQYYLDNPKERTEKIAECKNNPGEKALAANCVNATEAQRKAMFQGD  GMPKIR"
R388_heat.insert(-6, 0)
R388_heat.insert(-6, 0)
EEX_COMB = []
for i, j in zip(list(pMBUI4EEX), list(R388EEX)):
    EEX_COMB.append(i + '\n' + j)


xaxis_EEX = range(0, len(pMBUI_heat))
xaxis_EEX = [i + 0.5 for i in xaxis_EEX]

prev= ""

num=0

VirB5_R388 = [0]*len(R388VirB5_aligned)
VirB6_R388 = [0]*len(R388VirB6_aligned)
VirB5_pMBUI = [0]*len(pMBUIVirB5_aligned)
VirB6_pMBUI = [0]*len(pMBUIVirB6_aligned)

Vir_heat = [VirB5_pMBUI, VirB6_pMBUI, VirB5_R388, VirB6_R388]

for i in range(0, len(interactome["EEX"])):
    scores = [0,0,0,0]
    if interactome["EEX"][i] != prev:
        vir_chain = interactome["Vir"][i]
        print(vir_chain)
        lor, lop = [], []
        for j in fullist[num]:
            if vir_chain + ":" in j:
                pos = re.search(r' [0-9]+', j)
                amino = re.search(r'[A-Z]{3} ', j)
                lop.append(int(pos.group()))
                lor.append(str(amino.group()).strip())
        for w in Vir_list:
            for k in range(0,len(lor)): #Hay que tener en cuenta el alineamiento, que no corresponde con el valor de lor[k]
                try:
                    if resdict[lor[k]] == w[lop[k]-1]:
                        scores[Vir_list.index(w)] += 1
                except IndexError:
                    continue
        Vir_exact = scores.index(max(scores))
        for k in range(0,len(lop)):
            Vir_heat[Vir_exact][lop[k]] += 1
        num += 1
    prev = interactome["EEX"][i]

Vir_heat_aligned = [[], [], [], []]
for i in range(0, len(Vir_heat)):
    offset = 0
    for j in range(0, len(Vir_list_aligned[i])):
        if Vir_list_aligned[i][j] == Vir_list[i][j - offset]:
            Vir_heat_aligned[i].append(Vir_heat[i][j])
        else:
            Vir_heat_aligned[i].append(0)
            offset += 1

xaxis_VirB5 = range(0, len(Vir_heat_aligned[0]))
xaxis_VirB5 = [i + 0.5 for i in xaxis_VirB5]

xaxis_VirB6 = range(0, len(Vir_heat_aligned[1]))
xaxis_VirB6 = [i + 0.5 for i in xaxis_VirB6]

VirB5_COMB = []
for i, j in zip(list(Vir_list_aligned[0]), list(Vir_list_aligned[2])):
    VirB5_COMB.append(i + '\n' + j)

VirB6_COMB = []
for i, j in zip(list(Vir_list_aligned[1]), list(Vir_list_aligned[3])):
    VirB6_COMB.append(i + '\n' + j)

plt.xticks(xaxis_VirB5, VirB5_COMB)
plt.stairs(Vir_heat_aligned[2], alpha=0.5, label='R388 VirB5', fill = True)
plt.stairs(Vir_heat_aligned[0], alpha=0.5, label='pMBUI VirB5', fill = True)

plt.xlabel('Residue Position')
plt.ylabel('Frequency')
plt.title('Heatmap of Interactions')
plt.legend()
plt.show()

plt.xticks(xaxis_VirB6, VirB6_COMB)
plt.stairs(Vir_heat_aligned[3], alpha=0.5, label='R388 VirB6', fill = True)
plt.stairs(Vir_heat_aligned[1], alpha=0.5, label='pMBUI VirB6', fill = True)

plt.xlabel('Residue Position')
plt.ylabel('Frequency')
plt.title('Heatmap of Interactions')
plt.legend()
plt.show()



plt.xticks(xaxis_EEX, EEX_COMB)
plt.stairs(R388_heat, alpha=0.5, label='R388', fill = True)
plt.stairs(pMBUI_heat, alpha=0.5, label='pMBUI', fill = True)

plt.xlabel('Residue Position')
plt.ylabel('Frequency')
plt.title('Heatmap of Interactions')
plt.legend()
plt.show()

results.close()