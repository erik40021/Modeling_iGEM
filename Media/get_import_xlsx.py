from csv import writer
from linecache import getline
import xlsxwriter

def search(list, platform):
    for i in range(len(list)):
        if list[i] == platform:
            return True
    return False

YL=open("Data/iYLI647.xml","r")
lines=YL.readlines()
workbook = xlsxwriter.Workbook('Output/Transport_Reactions.xlsx')      #create .xlsx
worksheet = workbook.add_worksheet()
List=[]
cnt=0
worksheet.write(cnt,0,"reactions")
cnt+=1
readline=False
for line in lines: 
    
    if line.find("</reaction>")>0:
        readline=False
    if line.find("<p>SUBSYSTEM: Transport Extracellular</p>")>0:
        readline=True
    if(line.find("<speciesReference species=")and readline==True):
        part=line[line.find("<speciesReference species=")+27:]
        if(part.find("_e\"")>0):
            part=part[:part.find("stoichiometry=")-2]
            List.append(part)
        
Rxns=[]
for name in List:
    if search(Rxns,name) == False:
        Rxns.append(name)

for rx in Rxns:
    worksheet.write(cnt,0,rx)
    cnt+=1

workbook.close()


