from pylab import *
from DavidLib import load_xvg_Data, smooth


fidpi  = 100
fisize = (10,36)
fontsz = 10
m_lw   = 1.5
m_size = 4 


###matplotlib parameters
rcParams['font.sans-serif']='Arial'
rcParams['legend.fontsize']= fontsz-2
rcParams['font.size']= fontsz
rcParams['lines.markersize']= m_size

fidpi  = 1200
sadpi  = 1200
fisize = array([7.2,4.8])*1



Dists = [5,6,7,8,9,10,11,12]
Dirs  = [ "WT/WT_Knone/%02d_NS/"%i for i in Dists]

RMSF_1  = ["RMSF_O_%s_1.xvg"%s for s in ["A","B","C","D"] ] 
RMSF_2  = ["RMSF_O_%s_2.xvg"%s for s in ["A","B","C","D"] ]
RMSF_3  = ["RMSF_O_%s_3.xvg"%s for s in ["A","B","C","D"] ]
RMSF_4  = ["RMSF_O_%s_4.xvg"%s for s in ["A","B","C","D"] ]


T_sample = 60000

RMSF = load_xvg_Data("WT/WT_Knone/RMSF_Ebene.dat")*10
Proj_all = load_xvg_Data("WT/WT_Knone/Projections.dat")[:,0]






skip=1
Proj_N_05 = load_xvg_Data("WT/WT_Knone/05_NS/Filter_proj.xvg")[::skip*2]
Proj_F_05 = load_xvg_Data("WT/WT_Kfree/05_NS/Filter_proj.xvg")[::skip*8]
Proj_R_05 = load_xvg_Data("WT/WT_Kposre/05_NS/Filter_proj.xvg")[::skip*14]
Proj_N_10 = load_xvg_Data("WT/WT_Knone/10_NS/Filter_proj.xvg")[::skip]


Ichb = array([
	    (load_xvg_Data("WT/WT_Knone/ichb/NS_dist.xvg")[-54:,1]).mean()*10,
	    (load_xvg_Data("WT/WT_Knone/ichb/Filter_proj.xvg")[-54:,1]).mean()*10
	    ])
	    

D = load_xvg_Data("/home/dkoepfe/Desktop/Projects/Herg/Herg_Mutant_G628_S631C/Sim_Mutant_NoK/NS_dist.xvg")
Mut_dist = mean(D[ (D[:,0] > 40000) * (D[:,0] < 50000), 1 ])*10
D = load_xvg_Data("/home/dkoepfe/Desktop/Projects/Herg/Herg_Mutant_G628_S631C/Sim_Mutant_NoK/Filter_proj.xvg")
Mut_proj = mean(D[ (D[:,0] > 40000) * (D[:,0] < 50000), 1 ])*10

Mut = array([Mut_dist, Mut_proj])


figure(facecolor='w', figsize=fisize, dpi=fidpi)
subplot(2,2,2)
for c, R in zip(['ro-', 'go-', 'bo-', 'yo-'], RMSF.T):
  xlim(4.95,12.05)
  ylim(0.65, 1.7)
  yticks([0.8, 1.0, 1.2, 1.4, 1.6])
  plot(Dists, R[:len(Dists)], c, lw=m_lw, alpha =0.8)

plot(Dists, RMSF.mean(1), 'k--', lw=m_lw+2)

xlabel('Ser620 - Asn629 C  distance [A]')
ylabel('RMSF [A]')


subplot(2,2,3)
axhline(4.612, color = 'r', linewidth=m_lw, alpha = 0.8)
text(37.5, 4.8,'1K4D SF', horizontalalignment='center', verticalalignment='bottom', color = 'k')
axhline(-4.612, color = 'r', linewidth=m_lw, alpha = 0.8)
text(37.5, -4.9,'1K4C SF', horizontalalignment='center', verticaalignment='top', color = 'k')
xlim(-0.3,75); ylim(-6.5,6.5); yticks([-6, -4, -2, 0, 2, 4, 6])


sf = 30; nalp = 0.3; slw= m_lw; nlw=0.5; 

#c= 'm'
#plot(Proj_F_05[:, 0]/1000, Proj_F_05[:, 1]*10, c, alpha=nalp)
#plot(Proj_F_05[:, 0]/1000, smooth(Proj_F_05[:, 1]*10, sf), c, lw=slw, alpha=1, ls= "-")
#text(69, -1.4, '5A', color = c)

c= 'm'
plot(Proj_R_05[:, 0]/1000, Proj_R_05[:,1]*10, c, alpha=nalp, lw = 1)
plot(Proj_R_05[:, 0]/1000, smooth(Proj_R_05[:, 1]*10, sf*2), c, lw=slw, alpha=1, ls= ":")
text(69, -0.7, '5A', color = c)

c = 'b'
plot(Proj_N_05[:, 0]/1000, Proj_N_05[:, 1]*10, c, alpha=nalp, lw = 1)
plot(Proj_N_05[:, 0]/1000, smooth(Proj_N_05[:, 1]*10, sf), c, lw=slw)
text(60, 5, '5A', color = c)


c= 'g'
plot(Proj_N_10[:, 0]/1000, Proj_N_10[:,1]*10, c, alpha=nalp, lw = 1)
plot(Proj_N_10[:, 0]/1000, smooth(Proj_N_10[:, 1]*10, sf), c, lw=slw)
text(60, -5.7, '10A', color = c)

text(45, 0.2, "Ser620 - Asn629\nC  distance")


plot([],[],'k:', label='K+ in SF', lw=2)
plot([],[],'k-', label='without K+', lw=2)
legend(loc='upper left', shadow='True', fancybox='True')
ylabel('Pos. on interpolation vector [A]')
xlabel('Time [ns]')


subplot(2,2,4)
axhline(4.612, color = 'r', linewidth=m_lw, alpha = 0.8)
text(8.5, 4.8,'1K4D SF', horizontalalignment='center', verticalalignment='bottom', color = 'k')
axhline(-4.612, color = 'r', linewidth=m_lw, alpha = 0.8)
text(8.5, -4.9,'1K4C SF', horizontalalignment='center', verticalalignment='top', color = 'k')
ylim(-6.5,6.5); yticks([-6, -4, -2, 0, 2, 4, 6])
xlim(4.95,12.05)
plot(Dists, Proj_all*10, 'bo-', lw=m_lw)




plot(Ichb[0], Ichb[1], 'go', ms=m_size+2)
annotate('Inter-chain H-Bond', xy=(Ichb[0], Ichb[1]+.2), xycoords='data', xytext=(5, 35), textcoords='offset points', color='g',
	  arrowprops=dict(arrowstyle="simple", color='g',
                          fc="g", alpha=0.6, ec="none",
                          connectionstyle="arc3,rad=0.5"
                         )
        )

plot(Mut[0], Mut[1], 'ro', ms=m_size+2)
annotate('G628C-S631C\nMutant', xy=(Mut[0]+.1, Mut[1]+.2), xycoords='data', xytext=(10, 10), textcoords='offset points', color='r',
	  arrowprops=dict(arrowstyle="simple", color='r',
                          fc="r", alpha=0.6, ec="none",
                          connectionstyle="arc3,rad=0.3"
                         )
        )

xlabel('Ser620 - Asn629 C  distance [A]')
ylabel('Pos. on interpolation vector [A]')

subplots_adjust(left=0.06, bottom=0.1, top=1-0.01, right=1-0.01, wspace=0.3, hspace=0.3)
#subplots_adjust(hspace=0.1,wspace=0.1)
savefig('Figure1_mat.png', dpi=sadpi)

show()






  
  
      
