# coding: UTF-8
import sys
l1111l1 = sys.version_info [0] == 2
l1ll111l1 = 2048
l111l1 = 7
def l1ll1l1 (l1l111l1):
    global l1ll11l1
    l1l1l1 = ord (l1l111l1 [-1])
    l1llll1l1 = l1l111l1 [:-1]
    l11l1 = l1l1l1 % len (l1llll1l1)
    l1l11l1 = l1llll1l1 [:l11l1] + l1llll1l1 [l11l1:]
    if l1111l1:
        l11111l1 = unicode () .join ([unichr (ord (char) - l1ll111l1 - (l1l1l1l1 + l1l1l1) % l111l1) for l1l1l1l1, char in enumerate (l1l11l1)])
    else:
        l11111l1 = str () .join ([chr (ord (char) - l1ll111l1 - (l1l1l1l1 + l1l1l1) % l111l1) for l1l1l1l1, char in enumerate (l1l11l1)])
    return eval (l11111l1)
import matplotlib.pyplot as plt
def l1lll1l1(l1lll11l1):
    if isinstance(l1lll11l1,tuple):
        return max(l1lll1l1(l1lll11l1[0]),l1lll1l1(l1lll11l1[1]))+1
    else:
        return 0
def ll1l1(l1lll11l1,l111l1l1,l1ll1l1l1,l11l11l1=None):
    if not isinstance(l1lll11l1,tuple):
        x = l111l1l1
        y = l11l11l1[l1lll11l1] if l11l11l1 else l1ll1l1l1
        plt.plot([x,x],[0,y],l1ll1l1 (u"ࠫࡰ࠭ࠀ"))
        plt.text(x,0,l1lll11l1)
        return x,y,1
    else:
        xl,l11ll1l1,nl = ll1l1(l1lll11l1[0],l111l1l1,l1ll1l1l1-1,l11l11l1)
        xr,yr,nr = ll1l1(l1lll11l1[1],l111l1l1+nl,l1ll1l1l1-1,l11l11l1)
        y = l11l11l1[l1lll11l1] if l11l11l1 else l1ll1l1l1
        plt.plot([xl,xl],[l11ll1l1,y],l1ll1l1 (u"ࠬࡱࠧࠁ"))
        plt.plot([xr,xr],[yr,y],l1ll1l1 (u"࠭࡫ࠨࠂ"))
        plt.plot([xl,xr],[y,y],l1ll1l1 (u"ࠧ࡬ࠩࠃ"))
        return (xl+xr)/2,y,nl+nr
def l11l1l1(l1lll11l1,l11l11l1):
    ll1l1(l1lll11l1,0,l1lll1l1(l1lll11l1),l11l11l1)
    axes = plt.gca()
    lim = axes.get_ylim()
    delta=(lim[1]-lim[0])*0.05
    axes.set_ylim(lim[0]-delta,lim[1]+delta)
    lim = axes.get_xlim()
    delta=(lim[1]-lim[0])*0.05
    axes.set_xlim(lim[0]-delta,lim[1]+delta)
    plt.show()
def showtree(node,height=None):
	l11l1l1(node,height)