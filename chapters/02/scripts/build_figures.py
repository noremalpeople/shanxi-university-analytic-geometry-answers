from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(__file__).resolve().parent.parent/'figures'
out.mkdir(exist_ok=True)

def view(X,Y,Z,name,title,lim=(-3.5,3.5)):
    fig=plt.figure(figsize=(5.5,4.2))
    ax=fig.add_subplot(111,projection='3d')
    ax.plot_surface(X,Y,Z,alpha=0.64,linewidth=0.0,antialiased=True)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_zlim(*lim)
    ax.set_title(title,fontsize=12)
    ax.set_box_aspect((1,1,0.85))
    ax.view_init(elev=20,azim=-50)
    fig.tight_layout()
    fig.savefig(out/name,dpi=185,bbox_inches='tight')
    plt.close(fig)

u=np.linspace(0,2*np.pi,70); v=np.linspace(0,2*np.pi,32)
U,V=np.meshgrid(u,v)
R,r=2,0.55
X=(R+r*np.cos(V))*np.cos(U);Y=(R+r*np.cos(V))*np.sin(U);Z=r*np.sin(V)
view(X,Y,Z,'torus.png','Torus: rotating an offset circle')

u=np.linspace(0,2*np.pi,70); v=np.linspace(-2.0,2.0,44)
U,V=np.meshgrid(u,v)
X=np.sqrt(1+V**2)*np.cos(U);Y=np.sqrt(1+V**2)*np.sin(U);Z=V
view(X,Y,Z,'one-sheet-hyperboloid.png','One-sheet hyperboloid')

u=np.linspace(-1.9,1.9,60);v=np.linspace(-1.9,1.9,60)
X,Y=np.meshgrid(u,v);Z=0.42*(X**2-Y**2)
view(X,Y,Z,'hyperbolic-paraboloid.png','Hyperbolic paraboloid (saddle)')

u=np.linspace(0,2*np.pi,70);v=np.linspace(-2.4,2.4,70)
U,V=np.meshgrid(u,v)
X=0.72*np.abs(V)*np.cos(U);Y=0.72*np.abs(V)*np.sin(U);Z=V
view(X,Y,Z,'double-cone.png','Double cone and axis')

u=np.linspace(0,2*np.pi,70);v=np.linspace(0,np.pi,40)
U,V=np.meshgrid(u,v)
X=2.0*np.sin(V)*np.cos(U);Y=1.5*np.sin(V)*np.sin(U);Z=1.4*np.cos(V)
view(X,Y,Z,'ellipsoid.png','Ellipsoid (schematic)',lim=(-2.8,2.8))
print('figures:',[p.name for p in out.glob('*.png')])
