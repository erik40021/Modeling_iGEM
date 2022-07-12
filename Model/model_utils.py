# Useful functions to interact with a COBRApy model

# Don't forget: for almost every thing we want to do, there might already 
# be a perfectly fitting function!

from inspect import stack
import math
from matplotlib import pyplot as plt
import numpy as np
import xlsxwriter


def plot_fluxes(res, labels, colors, stack_indices, bar_width=0.1, size=(10,7), save_as=None):
    ''' 
    Plots selection of fluxes of solution vector as multiple bar chart.
    :param list res: 2-dimensional list of fluxes after each manipulation to the model. Each sublist
        must contain the same number of selected fluxes (in the same order) and the respective change-info at the end.
    :param list labels: list of labels of the selected fluxes.
    :param list colors: list of colors of the selected fluxes.
    :param stack_indices: list of 2 tuples. 1st tuple: Indices of fluxes in res to be stacked. 2nd tuple: Indices of fluxes on which
        to stack them. Example: [(5,6), (2,2)] means that on flux 2 there shall be stacked fluxes at indices 5 and 6 in res.
    :param float bar_width: width of all bars in the chart.
    :param tuple size: size of the chart, x for width, y for height.
    :param str save_as: optional, name of the file when saving is desired.
    '''
    subtitles = [i[len(res[0])-1] for i in res]
    ys = []
    for i in range(0, len(res[0])-1):
        ys += [j[i] for j in res]
    for i in range(0, len(res)):
        ys[len(res)*2 + i] += ys[len(res)*5 + i]
        ys[len(res)*3 + i] += ys[len(res)*6 + i] + ys[len(res)*7 + i]
    max_y = max(ys)
    for i in range(0,len(subtitles)):
        if len(subtitles[i]) > 14:
            a = str.split(subtitles[i])
            subtitles[i] = ' '.join(a[:math.ceil(len(a)/2)]) + "\n" + ' '.join(a[math.ceil(len(a)/2):])
    x = np.arange(len(subtitles))  # the label locations
    width = bar_width # the width of the bars
    fig, ax = plt.subplots()
    offset = -2
    for i in range(0, len(res[0])-(1+len(stack_indices[0]))):
        ax.bar(x + (offset + i)*width, [j[i] for j in res], width, label=labels[i], color=colors[i])
        if i in stack_indices[1]:
            index_prev = i
            for k in range(stack_indices[1].index(i), stack_indices[1].index(i)+stack_indices[1].count(i)):
                ax.bar(x + (offset + i)*width, [j[stack_indices[0][k]] for j in res], width, bottom=[j[index_prev] for j in res],
                    label=labels[stack_indices[0][k]], color=colors[stack_indices[0][k]])
                index_prev = stack_indices[0][k]
        
     # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Flux [mmol/gdcw/h]')
    ax.set_title('Fluxes after Manipulations')
    if len(res) > 5:
        ax.set_xticks(x, subtitles, rotation=45)
    else:
        ax.set_xticks(x, subtitles)
    ax.set_ylim(top=max_y + 0.05*max_y)
    ax.legend()
    fig.tight_layout()
    fig.set_size_inches(size[0], size[1])
    if save_as is not None:
        plt.savefig("Output/" + save_as, dpi=500)
    plt.show()


def write_to_excel(res, file_name, labels=None):
    '''
    Writes given fluxes to Excel file.
    :param list res: 2-dimensional list of fluxes. Each sublist must contain the same number of selected fluxes 
        (in the same order) and the respective change-info.
    :param str file_name: name of the Excel file to save the fluxes.
    :param list labels: optional, list of labels of the selected fluxes.
    '''
    if labels is None:
        labels = ["biomass", "fpp_prod", "aps_gpp_apinene", "erg20_gpp", "npp_prod", "aps_npp_apinene", "aggpps2_gpp", "mfps144_gpp"]
    results = [labels] + res
    with xlsxwriter.Workbook('Output/' + file_name + '.xlsx') as workbook:
        worksheet = workbook.add_worksheet()
        for row_num, data in enumerate(results):
            worksheet.write_row(row_num, 0, data)


def summarize_properties(model, short=False):
    """ informs about model properties"""
    print(f'\nModel summary:\n{len(model.reactions)} reactions')
    print(f'{len(model.metabolites)} metabolites')
    print(f'{len(model.genes)} genes')

    if not short:
        # Iterate through the objects in the model (details)
        max_entries_per_object = 10
        print(f"\nReactions (total: {len(model.reactions)})")
        print("---------")
        print("First %s reactions:" % max_entries_per_object)
        for x in model.reactions[:max_entries_per_object]:
            print("%s : %s" % (x.id, x.reaction))

        print(f"\nMetabolites (total: {len(model.metabolites)})")
        print("-----------")
        print("First %s metabolites:" % max_entries_per_object)
        for x in model.metabolites[:max_entries_per_object]:
            print('%9s : %s' % (x.id, x.formula))

        print(f"\nGenes (total: {len(model.genes)})")
        print("-----")
        print("First %s genes:" % max_entries_per_object)
        for x in model.genes[:max_entries_per_object]:
            associated_ids = (i.id for i in x.reactions)
            print("%s is associated with reactions: %s" %
                (x.id, "{" + ", ".join(associated_ids) + "}"))