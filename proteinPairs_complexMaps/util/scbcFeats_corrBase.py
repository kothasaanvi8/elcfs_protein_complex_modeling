import pandas as pd
import pickle as p

savedData = p.load(open('/home/gwilkins/manuscript_fig1_whiteboard_update3.pkl', 'rb'))
for key, item in savedData.items():
    globals()[key] = item

del key
del item
del savedData
for frac in scbcFractionation_data_byFraction.keys():
    print('nonspecific: ' + frac)
    print('generating correlation matrix...')
    placeholder = \
        scbcFractionation_data_byFraction[frac].T.corr().rename(
            columns=geneIDmapping,
            index=geneIDmapping).rename_axis(None).rename_axis(None, axis=1).stack().reset_index()
    placeholder = placeholder.rename(columns={'level_0': 'idi', 'level_1': 'idii', 0: 'corr'})
    print('appending unique pair IDs for subsequent duplicate removal...')
    placeholder['pairsFrozen'] = placeholder.loc[:, ['idi', 'idii']].apply(frozenset, axis=1)
    print('removing duplicates...')
    placeholder.drop_duplicates(subset=['pairsFrozen'], inplace=True)
    print('saving file...')
    p.dump(placeholder, open('/home/gwilkins/nonspecific_' + frac + '.pkl', 'wb'))

    for cell, cellSpecific_dict in scbcFractionation_data_byFraction_cellSpecific.items():
        print('specific: ' + cell)
        print('generating correlation matrix...')
        placeholder = \
            cellSpecific_dict[frac].T.corr().rename(
                columns=geneIDmapping,
                index=geneIDmapping).rename_axis(None).rename_axis(None, axis=1).stack().reset_index()
        placeholder = placeholder.rename(columns={'level_0': 'idi', 'level_1': 'idii', 0: 'corr'})
        print('appending unique pair IDs for subsequent duplicate removal...')
        placeholder['pairsFrozen'] = placeholder.loc[:, ['idi', 'idii']].apply(frozenset, axis=1)
        print('removing duplicates...')
        placeholder.drop_duplicates(subset=['pairsFrozen'], inplace=True)
        print('saving file...')
        p.dump(placeholder, open('/home/gwilkins/nonspecific_' + cell + '_' + frac + '.pkl', 'wb'))
