#!/usr/bin/env ipython
## This python script along with network_props.csv and polarization.csv will generate Figure 1 of the Nature Comms Paper
import graph_tool.all as gt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{mathptmx}']
matplotlib.rcParams['text.usetex'] = True

network_class = ['ws']
ix = pd.IndexSlice
colors = ['darkslateblue','darkcyan','coral']
#Toggle hatch
hatching = False
## Pick which cascade sizes are consider, valid choises are 0.1,0.2,...,0.9
cascades = list(map(str,list(np.round(np.linspace(0.1,0.9,9),1))))

network_props=pd.read_csv('network_props.csv')
network_props.set_index(['ID','network','p'],inplace=True)

polarization=pd.read_csv('polarization.csv')
polarization.set_index(['ID','network','p','th','seed'],inplace=True)

mpol = polarization.groupby(['p','th','network']).mean()
#Get mean polarization on network from seeds
polarization_mean = polarization.groupby(level=[0,1,2,3]).mean()

print('plotting')
save = True

hatch_regions = {'simple':{'x':[-10,polarization.index.get_level_values(3).unique()[3]],'y1':[-10,-10],'y2':[10,10],'hatch':'////','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10},
'trans':{'x':[polarization.index.get_level_values(3).unique()[3],polarization.index.get_level_values(3).unique()[8]],'y1':[-10,-10],'y2':[10,10],'hatch':'_','facecolor':'w','alpha':0,'edgecolor':'black','linewidth':1.0,'zorder':-10 },
'complex':{'x':[polarization.index.get_level_values(3).unique()[8],10],'y1':[-10, -10],'y2':[10,10],'hatch':'\\\\\\\\ ','facecolor':'w','alpha':0.2,'edgecolor':'black','linewidth':1.0,'zorder':-10}}

#For all values
pick_props = {'ws':[0,1,2,3,4,5,6,7,8,9,10]}
#For the 4 values used in the paper
pick_props = {'ws':[0,5,7,8]}

g = {'param1':[],'param2':[],'th':[]}
network_corr = pd.DataFrame(data=g)
network_corr.set_index(['param1','param2','th'],inplace=True)
methods = ['pearson','spearman']
flag_l_label=0
for network in network_class:
    probabilities = np.sort(mpol.loc[ix[:,:,network],:].index.get_level_values(0).unique())[pick_props[network]]
    for cas in cascades[:]:
        fig,axs = plt.subplots(figsize=(5.5*1.2,2.31*0.9),ncols=2,nrows=1,sharex=False,sharey=False)
        fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=None)
        axs[0].get_xaxis().get_major_formatter().set_scientific(False)
        ### Plotting loop
        axs[0].set_title(f'Cascade size {cas}')
        corr_fig_legend_label = [r'$C$',r'$\ell$',r'$R_g$']
        for idx,p in enumerate(probabilities[::-1]):
            C_label = str(network_props.loc[ix[:,network,p],:]['CC'].mean().round(2))
            if len(C_label) < 4:
                C_label =  C_label + '0'
            l_label = str(network_props.loc[ix[:,network,p],:]['Sp'].mean().round(1))
            if len(l_label) > 3:
                l_label = l_label[:3]
                flag_l_label=1
            print(len(l_label))

            r_label = str((network_props.loc[ix[:,network,p],:]['Rg'].mean()/1000).round(1))
            if len(r_label) > 3:
                r_label = r_label[:2]
            r_label = r_label + r'\! \times \! 10^{3}'

            if flag_l_label:
                label_string =r'${}  |  {}  |  {} $'.format(C_label,l_label,r_label)
            else:
                label_string =r'${}  |  {} \, |  {} $'.format(C_label,l_label,r_label)

            pol_fig_legend_label = label_string
            #Make lines distingushable
            if idx < 4:
               axs[0].plot(mpol.loc[ix[p,:,network],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network],f'{cas}'],ls='-',label=pol_fig_legend_label)
            elif idx < 8:
               axs[0].plot(mpol.loc[ix[p,:,network],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network],f'{cas}'],ls='.-',label=pol_fig_legend_label)
            elif idx < 12:
                axs[0].plot(mpol.loc[ix[p,:,network],f'{cas}'].index.get_level_values(1),mpol.loc[ix[p,:,network],f'{cas}'],ls='--',label=pol_fig_legend_label)

        # Set scale stuff
        axs[0].set_yscale('log')
        axs[0].set_ylabel(r'Polarization Speed $(v)$',labelpad=2.5)
        axs[0].set_xlabel(r'Threshold $( \theta )$',labelpad=2,math_fontfamily='cm')
        axs[0].set_ylim([1.95*10**-4,0.999999])
        axs[0].set_xlim([-0.02,0.54])
        ## Legend and title
        legend0 = axs[0].legend(title=r' $  C \;\, | \;\, \ell  \;\:  |\;\;\, R_{g} $', framealpha=1, facecolor='white',loc=[0.535,0.3],edgecolor='w',borderpad=0.2,markerscale=0.8,handlelength=1.4,handletextpad=0.4,fontsize=7)
        legend0.get_title().set_position((1.5,0))
        legend0.get_title().set_fontsize('7')
        for th in polarization.index.get_level_values(3).unique():
            (network_corr.loc[('CC','p',th),'speraman'],network_corr.loc[('Sp','p',th),'speraman'],network_corr.loc[('Rg','p',th),'speraman'])   = network_props.corrwith((polarization_mean.loc[ix[:,network,:,th],f'{cas}'].dropna().groupby('ID').mean()),method=methods[1])


        for idx,param in enumerate(network_corr.index.get_level_values(0).unique()):
            axs[1].plot(network_corr.loc[ix[param,:,:],:].index.get_level_values(2),network_corr.loc[ix[param,:,:],:],label=corr_fig_legend_label[idx],marker=['<','^','>'][idx%4],zorder=[3,2,1,4][idx%4],c=colors[idx%4],linewidth=2,markersize=[5,5,5][idx%4],ls=[':','--','-'][idx%4],markeredgewidth=[1,1,1][idx%4],markerfacecolor='none')

        ## set scale stuff
        axs[1].set_ylabel(r'Correlation $r_s(v,\chi)$',labelpad=2.5)
        axs[1].set_xlabel(r'Threshold $( \theta )$',labelpad=2)
        axs[1].set_ylim([-1.05,1.475])
        axs[1].set_xlim([-0.02,0.51])

        axs[0].text(-0.02,1.03,r'\textbf{(a)}',transform=axs[0].transAxes,fontsize=8)
        axs[1].text(-0.02,1.03,r'\textbf{(b)}',transform=axs[1].transAxes,fontsize=8)
        legend1 = axs[1].legend(title=r'$\chi$',loc=[0.7,0.2],borderpad=0.2,markerscale=1,handlelength=2,handletextpad=0.4,fontsize=7)

        legend1.get_title().set_position((0,0))
        legend1.get_title().set_fontsize('7')

        box_props = dict(alpha=1,facecolor='w',linewidth=0,zorder=100000,boxstyle='round',pad=0.5)
        ##Set background color to transitions
        if hatching:
            for x in hatch_regions:
                for ax in axs:
                    ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5])
                    ax.fill_between(**hatch_regions[x])
                    ax.tick_params(axis='x',labelsize=7)
                    ax.tick_params(axis='y',labelsize=7)

        for ax in axs:
            ax.text(0.055,0.9,r'{\fontfamily{phv}\selectfont  \textbf{Simple}}',transform=ax.transAxes,bbox=box_props,fontsize=7,fontdict={'family':'sans-serif'})
            ax.text(0.69,0.9,r'{\fontfamily{phv}\selectfont   \textbf{Complex}}',transform=ax.transAxes,bbox=box_props,fontsize=7,fontdict={'family':'sans-serif'})
        axs[0].grid(False)
        axs[1].grid(False)
        fig.show()
        fig.savefig(f'figures/fig1/fig1_{cas}.pdf')
