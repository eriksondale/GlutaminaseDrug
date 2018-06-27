from sys import argv

if(len(argv) != 3):
    print("Need input and output file....")
else:
    with open(argv[1],"r") as RxnFile:
        for line in RxnFile:
            line = line.strip("\n")
            line = line.split("\t")
            with open(argv[2],"a") as newRxnFile:
                rxn = line[1]
                newRxnFile.write(line[0] + "\t" + rxn[rxn.find('>>')+2:] + ">>" + rxn[0:rxn.find('>>')]+"\n")
RxnFile.close()
newRxnFile.close()
