from matplotlib import pyplot as plt
import numpy as np

def simString(
    t_max = 3, pluck_height = .02, pluck_pos = .5, 
    # sr = 44100, n_string_markers = 256, 
    sr = 44100, n_string_markers = 16, 
    spring_stiff = 1, string_density = 1, 
    wood_mass = 100,
):
    time_step = 1 / sr
    marker_mass = string_density / n_string_markers
    n_markers = n_string_markers + 2
    markers = np.zeros((n_markers, 2))
    markers_velocity = np.zeros((n_string_markers, 2))
    markers[:, 0] = np.linspace(0, 1, n_markers)
    pluc_i = round(n_markers * pluck_pos)
    markers[:pluc_i, 1] = np.linspace(
        0, pluck_height, pluc_i, 
    )
    markers[pluc_i:, 1] = np.linspace(
        pluck_height, 0, n_markers - pluc_i, 
    )
    wood_y = 0
    wood_velocity = 0
    clock_max = t_max * sr
    audio = np.zeros((clock_max, ))
    render(markers, n_markers, 0)
    for clock in range(clock_max):
        t = clock * time_step
        markers[0,  0] = 0
        markers[-1, 0] = 0
        markers[0,  1] = wood_y
        markers[-1, 1] = wood_y
        displace = markers[1:] - markers[:-1]
        wood_force = (displace[0, 1] - displace[-1, 1]) * spring_stiff
        wood_velocity += wood_force * (time_step / wood_mass)
        wood_y += wood_velocity * time_step
        forces = (displace[1:] - displace[:-1]) * spring_stiff
        markers_velocity += forces * (time_step / marker_mass)
        markers[1:-1] += markers_velocity * time_step
        audio[clock] = sampleAudio(wood_y)
        render(markers, n_markers, t)
    return audio

TSPF = 10   # time steps prt frame
tspf_progress = 0
def render(markers, n_markers, t, force = False):
    global tspf_progress
    if not force:
        tspf_progress -= 1
        if tspf_progress <= 0:
            tspf_progress = TSPF
        else:
            return
    plt.clf()
    plt.plot(markers[:, 0], markers[:, 1], '.')
    plt.axis('equal')
    plt.axis([0, 1, -.3, .3])
    plt.title(f't = {t}')
    plt.draw()
    plt.pause(.15)

def sampleAudio(wood_y):
    return wood_y

simString()
