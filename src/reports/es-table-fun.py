def display_table(t):
    display(Markdown(
        '<div class="datatable">' +
        tabulate(t, headers='keys') +
        "\n\n</div>"
    ))
def table_ds(ds, fdr, le_prop=0.0):
    df = ds.to_dataframe()
    df.reset_index(level=0, inplace=True)
    display_table(df.loc[(df['fdr'] < fdr) & (df['le_prop'] > le_prop)])
