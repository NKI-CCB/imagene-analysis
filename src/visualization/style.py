import matplotlib


def set_style(p=matplotlib.rcParams):
    p['font.family'] = 'sans-serif'
    p['font.sans-serif'] = ['Arial']
    p['font.weight'] = 'normal'
    p['figure.dpi'] = 300
    p['figure.facecolor'] = 'none'
    p['axes.labelsize'] = 11
    p['xtick.labelsize'] = 9
    p['ytick.labelsize'] = 9
