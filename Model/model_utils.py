# Useful functions to interact with a COBRApy model


# Don't forget: for almost every thing we want to do, there might already 
# be a perfectly fit function!



# (Definitely not complete) list of useful pre-installed functions:

# model.reaction.get_by_id()
# model.metabolite.get_by_id()



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