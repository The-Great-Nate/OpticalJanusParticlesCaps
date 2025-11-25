# Source - https://stackoverflow.com/a
# Posted by Joe Kington, modified by community. See post 'Timeline' for change history
# Retrieved 2025-11-10, License - CC BY-SA 3.0

import matplotlib.pyplot as plt
from matplotlib.patches import Wedge

def main():
    fig, ax = plt.subplots()
    dual_half_circle((0.5, 0.5), radius=0.3, angle=90, ax=ax)
    ax.axis('equal')
    plt.show()

def dual_half_circle(center, radius, angle=0, ax=None, colors=('w','k'),
                     **kwargs):
    """
    Add two half circles to the axes *ax* (or the current axes) with the 
    specified facecolors *colors* rotated at *angle* (in degrees).
    """
    if ax is None:
        ax = plt.gca()
    theta1, theta2 = angle, angle + 180
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[0], **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], **kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return [w1, w2]

main()
