import xlsxwriter

def medium_objectivevalue_xlsx(Food,Solutions):
    workbook = xlsxwriter.Workbook('Output/YL_growth_media.xlsx')      #create .xlsx
    worksheet = workbook.add_worksheet()
    worksheet.write(0,0,"This C-source is set to 1000")         #head of table
    worksheet.write(0,1,"objective value")
    row_num=0
    line=0
    for data in Solutions:
        if data>0:                                               #only print if solution is feasible
            worksheet.write(line+1,0, Food[row_num])             #write name of C-Source in first column
            worksheet.write(line+1,1, data)                      #write objective value for this c-source
            line+=1
        row_num+=1
    workbook.close()
    print("YL_growth_media.xlsx was written")
