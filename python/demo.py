import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import subprocess as sp

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, data):
        self.counter = 0
        self.data = data
        self.num_figs = len(data)
        self.stream = self.data_stream()
        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots()
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=10, 
                                          init_func=self.setup_plot, blit=True)

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        sim = np.array(self.data[0]).reshape((-1, 2))
        x, y = sim[:, 0], sim[:, 1]
        self.scat = self.ax.scatter(x, y, s=100)
                                    #cmap="jet", edgecolor="k")
        self.ax.axis([-1, 5, -3, 3])
        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,

    def data_stream(self):
        """Generate a random walk (brownian motion). Data is scaled to produce
        a soft "flickering" effect."""
        while True:
            sim = np.array(self.data[self.counter]).reshape((-1, 2))
            #x, y = sim[:, 0], sim[:, 1]
            self.counter += 1
            self.counter %= self.num_figs
            yield sim

    def update(self, i):
        """Update the scatter plot."""
        data = next(self.stream)
        # Set x and y data...
        self.scat.set_offsets(data)
        # Set sizes...
        #self.scat.set_sizes(300 * abs(data[:, 2])**1.5 + 100)
        # Set colors..
        #self.scat.set_array(data[:, 3])

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,


if __name__ == '__main__':
    binary = '../build/main'
    output_lines = sp.check_output([binary]).decode().split('\n')
    data = [[float(x) for x in l.split()] for l in output_lines if len(l) > 0]
    a = AnimatedScatter(data)
    a.ani.save('demo.gif', fps=50, dpi=80)
    plt.show()