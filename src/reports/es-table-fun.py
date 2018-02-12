def display_table(t):
    display(Markdown(
        '<div class="datatable">' +
        tabulate(t, headers='keys') +
        "\n\n</div>"
    ))
def table_ds(ds, fdr, le_prop=0.0, nes=0.0):
    df = ds.to_dataframe()
    df.reset_index(level=0, inplace=True)
    df = df.loc[(df['fdr'] < fdr) &
                (df['le_prop'] > le_prop) &
                (np.abs(df['nes']) > nes) ]
    display_table(df)
