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
workbook = xlsxwriter.Workbook('Output/EX_Reactions.xlsx')      #create .xlsx
worksheet = workbook.add_worksheet()
List=[]
cnt=0
worksheet.write(cnt,0,"reactions")
cnt+=1
for line in lines:
    if line.find("EX_")>0:
        part=line[line.find("EX_"):]
        List.append(part[:part.find("e_RPAREN_")+9])
        
Rxns=[]
for name in List:
    if search(Rxns,name) == False:
        Rxns.append(name)

for rx in Rxns:
    worksheet.write(cnt,0,rx)
    cnt+=1

workbook.close()


