# /usr/bin/env python

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt


# this will come in handy later
def distance(coord_1, coord_2):
    return np.sqrt((coord_1[0] - coord_2[0]) ** 2 + (coord_1[1] - coord_2[1]) ** 2)


class individual:
    # each individual is instantiated with a location,
    # a status and a set of boundaries
    def __init__(self, coords, status, bounds, speed):
        self.coord = coords
        self.status = status
        self.bounds = bounds
        self.speed = speed
        self.time = 0
        self.removed = 0

    # in this example I imagine each indivual moving in a given direction
    # which we can change at any point in time
    def set_direction(self):
        angle = np.random.uniform(0, 2 * np.pi)
        x_direction = np.cos(angle)
        y_direction = np.sin(angle)

        self.x_direction = x_direction
        self.y_direction = y_direction

    # Every turn, the individual will updat their location
    def update_location(self):
        speed = self.speed
        bounds = self.bounds

        distance = speed * np.random.random()
        # What to do if they run up against the boundary of the region?
        if self.coord[0] + distance * self.x_direction < 0:
            new_x = -distance * self.x_direction - self.coord[0]
            self.x_direction = -self.x_direction
        elif self.coord[0] + distance * self.x_direction > bounds[0]:
            new_x = 2 * bounds[0] - distance * self.x_direction - self.coord[0]
            self.x_direction = -self.x_direction
        else:
            new_x = self.coord[0] + distance * self.x_direction

        if self.coord[1] + distance * self.y_direction < 0:
            new_y = -distance * self.y_direction - self.coord[1]
            self.y_direction = -self.y_direction
        elif self.coord[1] + distance * self.y_direction > bounds[1]:
            new_y = 2 * bounds[1] - distance * self.y_direction - self.coord[1]
            self.y_direction = -self.y_direction
        else:
            new_y = self.coord[1] + distance * self.y_direction

        self.coord = (new_x, new_y)

    # each individual has a certain chance to be infected by others
    # within a certain radius
    def transmission(self, others, radius, probability):
        # of course if they're already recoved it's moot
        if self.removed == 1:
            pass
        # and if they're already infected, they won't get infected again
        elif self.status == 1:
            # but we can track how long they've been infected
            self.time += 1
            # And put them into the 'removed' category
            if self.time == 175:
                self.status = 2
                self.removed = 1
        else:
            # but if they're susceptible,
            neighbor_count = 0
            for other in others:
                # then every infected person nearby increases the probability of transmission
                if other.status == 1:
                    if distance(self.coord, other.coord) < radius:
                        neighbor_count += 1
            if sum([1 for x in np.random.random(neighbor_count) if x > probability]) > 0:
                self.status = 1


# setting the global variables up front
bounds = (10, 10)
radius = .2
probability = .2

everyone = []
# five infected individuals, at random coordinates
for n in range(5):
    everyone.append(individual((10 * np.random.random(), 10 * np.random.random()),
                               1, bounds, np.random.random()))
# 995 uninfected individuals
for n in range(995):
    everyone.append(individual((10 * np.random.random(), 10 * np.random.random()),
                               0, bounds, np.random.random()))
for person in everyone:
    person.set_direction()

infected = [5]
removed = [0]
# For each turn we'll
for n in range(45):
    # see if anyone is infected
    for i in range(0, len(everyone)):
        everyone[i].transmission(everyone[:i] + everyone[i + 1:], radius, probability)
    # update everyone's location
    for person in everyone:
        person.update_location()
    # keep track of how many people there are in each category
    infected.append(sum([x.status for x in everyone]))
    removed.append(sum([x.removed for x in everyone]))
vulnerable = [1000 - removed[i] - infected[i] for i in range(0, len(infected))]

# creating our animated figure
fig = plt.figure(figsize=(6, 4))
plt.xlim(0, 45)
plt.ylim(0, 1000)
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax = plt.axes()
# our masrterlist of artist objects
ims = []
for n in range(1, len(infected) + 1):
    x = list(range(0, len(infected)))[:n]
    y = [infected[:n], vulnerable[:n], removed[:n]]
    # each artist object is a stackplot figure
    im = ax.stackplot(x, y, labels=['infected', 'vulnerable', 'removed'],
                      colors=['r', 'b', 'grey'])
    ims.append(im)

# animate them using the artistanimation function
ani = animation.ArtistAnimation(fig, ims, interval=120, blit=True, repeat_delay=1000)
ani.save('animation.gif', writer='Pillow')
plt.show()
