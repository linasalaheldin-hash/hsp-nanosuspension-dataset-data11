import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
NAVY='#1d3b5f'; GR='#5a8a5a'; RED='#a3261a'
fig,ax=plt.subplots(figsize=(10,8.6)); ax.set_xlim(0,100); ax.set_ylim(0,86); ax.axis('off')
def box(x,y,w,h,fc,ec,t1,t2,c,bold=True):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=1.2',fc=fc,ec=ec,lw=2))
    ax.text(x+w/2+4,y+h*0.64,t1,ha='center',va='center',fontsize=11,fontweight='bold',color=c)
    ax.text(x+w/2+4,y+h*0.30,t2,ha='center',va='center',fontsize=9,color=c)
def icon(x,y,s,c):
    ax.add_patch(plt.Circle((x,y),2.6,color=c)); ax.text(x,y,s,ha='center',va='center',color='white',fontsize=13,fontweight='bold')
def arrow(x,y0,y1,c=NAVY):
    ax.annotate('',xy=(x,y1),xytext=(x,y0),arrowprops=dict(arrowstyle='-|>',color=c,lw=2.5))
def side(y,t1,t2,fc,ec,c):
    ax.add_patch(FancyBboxPatch((57,y),40,8,boxstyle='round,pad=0,rounding_size=1.2',fc=fc,ec=ec,lw=2))
    ax.text(77,y+5.2,t1,ha='center',va='center',fontsize=10,fontweight='bold',color=c)
    ax.text(77,y+2.4,t2,ha='center',va='center',fontsize=9.5,color=c)
    ax.annotate('',xy=(57,y+4),xytext=(50,y+4),arrowprops=dict(arrowstyle='-|>',color=ec,lw=2.5))
H=9; X=1; W=49
ys=[76,61,46,31,16,1]
box(X,ys[0],W,H,'#eaf1f8',NAVY,'Literature search & formulation curation','Drug + stabilizer · traceable DOI · outcome',NAVY); ax.add_patch(plt.Circle((6,ys[0]+4.5),2.6,color=NAVY)); ax.add_patch(plt.Circle((5.6,ys[0]+4.9),1.0,fill=False,ec='white',lw=1.8)); ax.plot([6.4,7.4],[ys[0]+4.1,ys[0]+3.1],color='white',lw=2)
box(X,ys[1],W,H,'#e4eee2',GR,'Formulations identified','n = 62  (57 unique APIs)',GR); icon(6,ys[1]+4.5,'✓',GR)
box(X,ys[2],W,H,'#eff6ea','#7aa07a','Remaining after Stage 1','n = 57','#3a5a3a'); icon(6,ys[2]+4.5,'→','#7aa07a')
box(X,ys[3],W,H,'#eff6ea','#7aa07a','Remaining after Stage 2','n = 43','#3a5a3a'); icon(6,ys[3]+4.5,'→','#7aa07a')
box(X,ys[4],W,H,'#e4eee2',GR,'Curated Stable subset (Table 1)','n = 41  (39 unique APIs)',GR); icon(6,ys[4]+4.5,'★',GR)
box(X,ys[5],W,H,'#dcebd9',GR,'Final analysis subset','n = 54  (42 unique APIs)  ·  43 Stable, 11 Unstable',GR); icon(6,ys[5]+4.5,'★',GR)
for i in range(5): arrow(25.5,ys[i],ys[i+1]+H)
side(ys[1]+0.5,'Excluded: no polymer HSP data','for the stabilizing agent · n = 5','#fbebe8','#c0503c',RED)
side(ys[2]+0.5,'Excluded · n = 14','4 HSP conflict · 1 shared source · 9 not recorded','#fbebe8','#c0503c',RED)
side(ys[3]+0.5,'Excluded: relative, not absolute, outcome','vs. a comparator stabilizer · n = 2','#fbebe8','#c0503c',RED)
ax.add_patch(FancyBboxPatch((57,ys[5]+9),40,8,boxstyle='round,pad=0,rounding_size=1.2',fc='#e4eee2',ec=GR,lw=2))
ax.text(77,ys[5]+14.2,'Added: failure-specific search',ha='center',va='center',fontsize=10,fontweight='bold',color=GR)
ax.text(77,ys[5]+11.4,'+13 rows (11 Unstable, 2 Stable) · 3 new APIs',ha='center',va='center',fontsize=9.5,color=GR)
ax.plot([57,25.5],[ys[5]+12.5,ys[5]+12.5],color=GR,lw=2.5)
fig.savefig('fig2_new.png',dpi=300,bbox_inches='tight',facecolor='white')
