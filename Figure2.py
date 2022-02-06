import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
from matplotlib.colors import to_hex
import matplotlib
import sys
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{mathptmx}']
matplotlib.rcParams['text.usetex'] = True

ix=pd.IndexSlice

networks = ['smallworld']
nodes = 240
degrees = [4,6,8,10,12,14,16,18,20,22,24,26,28,30,32]

def linepointstyle(order, color, marker):
    """Return style parameters for a
    solid line of a given color and
    order with points denoted by a
    given marker.
    """
    return {'zorder': order, 'lw': 1.5, 'c': color, 'ls': '-', 'marker': marker,
            'markeredgecolor': color, 'markersize': 5, 'clip_on': False}

def pointstyle(order, color, marker):
    """Return style parameters for a
    points of a given color and order,
    denoted by a given marker.
    """
    return {'zorder': order, 'color': color, 'marker': marker,
            's': 30, 'clip_on': False}

colors = ['r','b','g','y']
markers = 's ^ o v'.split()
lp_main = [linepointstyle(10 - i, colors[i], markers[i]) for i in range(4)]


net_styles = [{'node_color': to_hex(p['c']),
               'node_shape': '.',
               'node_shape': p['marker'],
               'width': 0.005, #p['lw'],
               'edge_color': to_hex(p['c']),
               'node_size': 0.01} for p in lp_main]


ps_main = [pointstyle(10 - i, colors[i], markers[i]) for i in range(4)]

print('Loading data ..,')

if True:
    d = {'ID':[],'freq':[],'p':[]}
    network_gains = pd.DataFrame(data=d)
    network_gains.set_index(['ID','freq','p'],inplace=True)

    e = {'ID':[],'p':[]}
    network_props = pd.DataFrame(data=e)
    network_props.set_index(['ID','p'],inplace=True)

    insert = True
    landscape = True

    main_props = ['CC','SP','Rg']
    aux_props = ['rawCC','rawSP','rawRg','k','p','T']
    props = main_props + aux_props
    prop_label= {'CC':r'$\bar{C}$','SP':r'$\bar{\ell}$','Rg':r'$\bar{R}_g$'}
    for network in networks:
        for k in degrees:
            print(k)
            new_network_gains = pd.read_csv(f'networks/{nodes}/{k}/{network}_corr_gains.csv',sep='\t',index_col=[1])
            new_props = pd.read_csv(f'networks/{nodes}/{k}/{network}_props.csv',sep='\t',index_col=[0,1])

            new_props['p'] = new_props.index.get_level_values(1)

            new_props['k'] = k
            new_props['rawCC'] = new_props.CC
            new_props['rawSP'] = new_props.SP
            new_props['rawRg'] = new_props.Rg
            new_props.CC = new_props.CC/new_props.CC.max()
            new_props.loc[:,'T'] = new_props.loc[:,'T']/new_props.loc[:,'T'].max()
            new_props.Rg = new_props.Rg/new_props.Rg.min()
            new_props.SP = new_props.SP/new_props.SP.min()

            for f in new_network_gains.index.unique():
                new_network_gains.loc[f,'normH2'] = new_network_gains.loc[f].H2/new_network_gains.loc[f,'H2'].max()


            if 'k' in new_network_gains.columns:
                new_network_gains.loc[network_gains.isnull().k,'k']=k
            else:
                new_network_gains['k'] = k

            network_props = network_props.append(new_props.loc[:,props])

            new_network_gains = new_network_gains.reset_index()
            new_network_gains['freq'] = new_network_gains.freq.apply(lambda x: round(x,5))

            new_network_gains.set_index(['ID','freq'],inplace=True)

            network_gains = network_gains.append(new_network_gains)

    network_props.fillna(value=0,inplace = True)
    network_gains.fillna(value=0,inplace = True)

    print('Getting correlations ..,')
    freq16 = network_gains.index.get_level_values(1).unique()
    freqAll = network_gains.loc[network_gains.k==[x for x in degrees if x != 16][0]].index.get_level_values(1).unique()

    gg = pd.MultiIndex.from_tuples(list(zip(main_props*len(freqAll),['H2']*len(freqAll)*len(main_props),sorted(list(freqAll)*len(main_props)))))
    corr_index_16 = pd.MultiIndex.from_tuples(list(zip(main_props*len(freq16),['H2']*len(freq16)*len(main_props),sorted(list(freq16)*len(main_props)))))
    network_corr_16 = pd.DataFrame(index=corr_index_16)
    network_corr_lim = pd.DataFrame(index=gg)

    correlations = ['spearman']

    max_sp_lim = network_props.SP.max()
    min_sp_lim = network_props.SP.min()

    max_cc_lim = network_props.CC.max()
    min_cc_lim = network_props.CC.min()

    max_rg = network_props.Rg.max()
    min_rg = network_props.Rg.min()


    # # Good correlation for CC Rg
    min_sp_lim = 1.3
    max_sp_lim = 2.2

    network_props.query(f'{min_cc_lim} < CC < {max_cc_lim} and {min_sp_lim}  < SP < {max_sp_lim} and {min_rg} < Rg < {max_rg}').loc[:,('CC','SP','Rg')].corr(method='spearman').round(2)

    # query for  range of clutsering
    query_string = f'{min_cc_lim} < CC < {max_cc_lim} and {min_sp_lim} < SP < {max_sp_lim} and {min_rg} < Rg < {max_rg}'
    query_string16 = f'k == 16'

    lim_IDs = network_props.query(query_string).index.get_level_values(0)
    temp_network_ins = network_gains.loc[ix[lim_IDs,:]]
    corr_values_16 = []
    corr_values_lim = []
    for coef in correlations:
        ## Need to be sorted, because network_corr_... expects frequencies to be ordered
        for f in freq16.sort_values():
            corr_values_16 += list(network_props.query(query_string16).loc[:,main_props].corrwith(network_gains.loc[ix[:,f],'H2'],method=coef).values)

        for f in freqAll.sort_values():
            corr_values_lim += list(network_props.query(query_string).loc[:,main_props].corrwith(network_gains.loc[ix[:,f],'normH2'],method=coef).values)

network_corr_16.loc[corr_index_16,coef] = np.reshape(corr_values_16,(len(corr_values_16),1))
network_corr_lim.loc[gg,coef] = np.reshape(corr_values_lim,(len(corr_values_lim),1))
        
print('Plotting ..,')

hatch_regions = {'simple':{'x':[0.00005,network_corr_lim.loc['CC'].iloc[8].name[1]],'y1':[-1000,-1000],'y2':[10000,10000],'hatch':'////','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10},
'complex':{'x':[network_corr_lim.loc['CC'].iloc[9].name[1],100],'y1':[-2,-2],'y2':[10000,10000],'hatch':'\\\\\\\\ ','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10}}

my_new_colors = ['darkslateblue','darkcyan','coral']

cf = -3
sf = 0

specialp = [10,8,6,0]
specialp = [9,8,6,0]

fig,axs = plt.subplots(figsize=(5.5*1.2,2.1*1.7),ncols=2,nrows=2,sharex=False,sharey=False)
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.25)
## Getting correct plot setup
gs = axs[0][0].get_gridspec()
axs[0][0].remove()
axs[1][0].remove()
bigax = fig.add_subplot(gs[0:,0])

ax_corr_16 = axs[0][1]
ax_corr_lim = axs[1][1]

## Formating collective frequency response
bigax.set_xscale('log')
bigax.set_yscale('log')
bigax.set_xlabel(r'Frequency $(\omega)$',labelpad=2)
bigax.set_ylabel(r'Collective response $(H^2)$',labelpad=2.5)

specialk = 16
network_gains = network_gains.reset_index().set_index(['ID','freq','p'])
prop_special = network_props.loc[network_props.k==specialk].groupby(level=1).mean()
gain_special = network_gains.loc[network_gains.k==specialk].groupby(level=[2,1]).mean()
for p in gain_special.index.get_level_values(0).unique()[specialp]:
    C_label = str(prop_special.loc[p].rawCC.round(2))
    if len(C_label) < 4:
        C_label =  C_label + '0'
    l_label = str(prop_special.loc[p].rawSP.round(2))
    if len(l_label) > 3:
        l_label = l_label[:3]
        flag_l_label=1

    r_label = str((prop_special.loc[p].rawRg.mean()/1000).round(1))
    if len(r_label) > 3:
        r_label = r_label[:2]
    r_label = r_label + r'\! \times \! 10^{3}'

    if flag_l_label:
        label_string =r'${}  |  {}  |  {} $'.format(C_label,l_label,r_label)
    else:
        label_string =r'${}  |  {} \, |  {} $'.format(C_label,l_label,r_label)

    pol_fig_legend_label = label_string


    bigax.plot(gain_special.loc[p].H2.iloc[sf:cf],label=label_string)

leg1 = bigax.legend(title=r'$ C \;\, | \;\, \ell  \;\:  |\;\;\, R_{g} $', loc=[0.025,0.26],borderpad=0.2,markerscale=0.8,handlelength=0.9,handletextpad=0.4,fontsize=7)
leg1.get_title().set_position((3.55, 0))
leg1.get_title().set_fontsize('7')

# ,marker=['<','^','>'][idx%4]

for idx,metric in enumerate(main_props):
    ax_corr_16.plot(network_corr_16.loc[(metric,'H2')].iloc[sf:cf],label=prop_label[metric],zorder=[3,2,1,4][idx%4],c=my_new_colors[idx%4],linewidth=2,markersize=[5,5,5][idx%4],ls=[':','--','-'][idx%4],markeredgewidth=[1,1,1][idx%4],markerfacecolor='none')
    ax_corr_lim.plot(network_corr_lim.loc[(metric,'H2')].iloc[sf:cf],zorder=[3,2,1,4][idx%4],c=my_new_colors[idx%4],linewidth=2,markersize=[5,5,5][idx%4],ls=[':','--','-'][idx%4],markeredgewidth=[1,1,1][idx%4],markerfacecolor='none')
    

leg2=ax_corr_16.legend(title=r'$\bar{\chi}$',loc=[0.08,0.31],borderpad=0.2,markerscale=1,handlelength=2,handletextpad=0.4,fontsize=7)
leg2.get_title().set_position((0, 0))
leg2.get_title().set_fontsize('7')

insert_16 = fig.add_axes((0.74,0.62,0.155,0.1))

ps_main[3]['s'] = 0.5
ps_main[3]['zorder'] = 20
insert_16.scatter(network_props.query(query_string16).Rg,network_props.query(query_string16).CC,**ps_main[3])

insert_16.set_xlabel(r'$\bar{R}_g$',labelpad=-8)
insert_16.set_ylabel(r'$\bar{C}$',labelpad=-6)
insert_16.set_xlim([0.95,2.25])
insert_16.set_yticks([0,1])
insert_16.xaxis.set_label_coords(0.5,-0.1)
insert_16.set_ylim([0,1.02])
insert_16.set_xscale('log')
insert_16.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
insert_16.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([1,2]))
insert_16.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())

insert_lim = fig.add_axes((0.74,0.17,0.155,0.095))

insert_lim.scatter(network_props.query('k == 16 ').Rg,network_props.query(' k == 16').CC,**ps_main[3])

ps_main[2]['s'] = 1
insert_lim.scatter(network_props.query(' k != 16 ').Rg,network_props.query(' k != 16 ').CC,**ps_main[2])

rectangle = plt.Rectangle((network_props.query(query_string).Rg.min(),network_props.query(query_string).CC.min()),network_props.query(query_string).Rg.max()-network_props.query(query_string).Rg.min(),network_props.query(query_string).CC.max()-network_props.query(query_string).CC.min()+0.02,zorder=100,ec='red',fill=False)
insert_lim.add_patch(rectangle)

insert_lim.set_xlabel(r'$\bar{R}_g$',labelpad=-8)
insert_lim.set_ylabel(r'$\bar{C}$',labelpad=-5)
insert_lim.set_xscale('log')
insert_lim.set_xlim([0.95,10.35])
insert_lim.set_ylim([-0.03,1.03])
insert_lim.set_yticks([0,1])
insert_lim.xaxis.set_label_coords(0.55,-0.1)
insert_lim.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
insert_lim.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([1,2,5,10]))
insert_lim.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
## Setup axises
ax_corr_16.set_xscale('log')
ax_corr_16.set_ylabel(r'Correlation $r_s(\bar{H}^2(\omega),\bar{\chi})$',labelpad=2.5)

ax_corr_lim.set_xscale('log')
ax_corr_lim.set_xlabel(r'Frequency $( \omega )$',labelpad=2)

## Fancy shared y-label
ax_corr_16.yaxis.set_label_coords(-0.1,-0.1)

## Set axis limits
bigax.set_xlim([0.0001,2])
bigax.set_ylim([0.017,520])

ax_corr_16.set_ylim([-1.65,1.72])
ax_corr_16.set_xlim([0.0005,3])

ax_corr_lim.set_ylim([-1.85,1.73])
ax_corr_lim.set_xlim([0.0005,3])

## Indentifying letters
bigax.text(-0.02,1.03,r'\textbf{(a)}',transform=bigax.transAxes,fontsize=8)
ax_corr_16.text(-0.02,1.04,r'\textbf{(b)}',transform=ax_corr_16.transAxes,fontsize=8)
ax_corr_lim.text(-0.02,1.04,r'\textbf{(c)}',transform=ax_corr_lim.transAxes,fontsize=8)

ax_corr_lim.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([-1,0,1]))
ax_corr_16.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([-1,0,1]))

ax_corr_lim.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax_corr_16.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())

## Set hatch regions
for x in hatch_regions:
    bigax.fill_between(**hatch_regions[x]) 
    ax_corr_16.fill_between(**hatch_regions[x]) 
    ax_corr_lim.fill_between(**hatch_regions[x]) 
    hatch_regions[x].pop('hatch',None)
    bigax.fill_between(**hatch_regions[x]) 
    ax_corr_16.fill_between(**hatch_regions[x]) 
    ax_corr_lim.fill_between(**hatch_regions[x]) 

## text boxes
box_props = dict(alpha=1,facecolor='w',linewidth=0,zorder=100000,boxstyle='round',pad=0.4)

bigax.text(0.155,0.955,r'{\fontfamily{phv}\selectfont  \textbf{Simple}}',transform=bigax.transAxes,fontsize=7,bbox=box_props,fontdict={'family':'sans-serif'})
bigax.text(0.66,0.955,r'{\fontfamily{phv}\selectfont  \textbf{Complex}}',transform=bigax.transAxes,fontsize=7,bbox=box_props,fontdict={'family':'sans-serif'})


ax_corr_16.text(0.1,0.89,r'{\fontfamily{phv}\selectfont  \textbf{Simple}}',transform=ax_corr_16.transAxes,fontsize=7,bbox=box_props,fontdict={'family':'sans-serif'})
ax_corr_16.text(0.65,0.89,r'{\fontfamily{phv}\selectfont  \textbf{Complex}}',transform=ax_corr_16.transAxes,fontsize=7,bbox=box_props,fontdict={'family':'sans-serif'})

ax_corr_lim.text(0.1,0.89,r'{\fontfamily{phv}\selectfont  \textbf{Simple}}',transform=ax_corr_lim.transAxes,fontsize=7,bbox=box_props,fontdict={'family':'sans-serif'})
ax_corr_lim.text(0.65,0.89,r'{\fontfamily{phv}\selectfont  \textbf{Complex}}',transform=ax_corr_lim.transAxes,fontsize=7,bbox=box_props,fontdict={'family':'sans-serif'})

bigax.grid(False)
ax_corr_16.grid(False)
ax_corr_lim.grid(False)

bigax.tick_params(axis='x',labelsize=7)
ax_corr_16.tick_params(axis='x',labelsize=7)
ax_corr_lim.tick_params(axis='x',labelsize=7)
bigax.tick_params(axis='y',labelsize=7)
ax_corr_16.tick_params(axis='y',labelsize=7)
ax_corr_lim.tick_params(axis='y',labelsize=7)

insert_16.tick_params(axis='x',labelsize=7)
insert_lim.tick_params(axis='x',labelsize=7)
insert_16.tick_params(axis='y',labelsize=7)
insert_lim.tick_params(axis='y',labelsize=7)

random  = nx.connected_watts_strogatz_graph(20,8,1,seed=10)
lattice  = nx.connected_watts_strogatz_graph(20,8,0,seed=10)
network_scaling = 0.95
ran_ins = fig.add_axes((0.27,0.09,0.2*network_scaling,0.3*network_scaling))
lat_ins = fig.add_axes((0.3,0.45,0.2*1.1,0.3*1.1))

pos = []
step=2*np.pi/(np.max(lattice.nodes)+1)
for x in lattice.nodes:
     pos.append((np.sin(step*x),np.cos(step*x)))

poss = dict(zip(lattice.nodes,pos))
nx.draw_networkx_nodes(lattice,ax=lat_ins,pos=poss,node_size=15,node_color='#984ea3',node_shape='o',linewidths=0)
nx.draw_networkx_edges(lattice,ax=lat_ins,pos=poss,edge_color='#984ea3',width=0.5)

nx.draw_networkx_nodes(random,ax=ran_ins,pos=poss,node_size=15,node_color='#e41a1c',node_shape='o',linewidths=0)
nx.draw_networkx_edges(random,ax=ran_ins,pos=poss,edge_color='#e41a1c',width=0.5)

lat_ins.set_ylim([-1.5,1.5])
lat_ins.set_xlim([-1.9,1.9])
lat_ins.set_axis_off()
lat_ins.set_clip_on(False)

ran_ins.set_ylim([-1.5,1.5])
ran_ins.set_xlim([-1.9,1.9])
ran_ins.set_axis_off()
ran_ins.set_clip_on(False)

### White insert patch
ax_corr_16.fill_between(x=[0.0325,2.9],y1=[-1.35,-1.35],y2=[0.25,0.25],facecolor='w')
ax_corr_lim.fill_between(x=[0.021,2.9],y1=[-1.725,-1.725],y2=[-0.08,-0.08],facecolor='w',zorder=-10)


fig.show()
fig.savefig('figures/fig2/transition_in_LFC.pdf')
