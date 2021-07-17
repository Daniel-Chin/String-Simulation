'''
spring rest length
'''
from matplotlib import pyplot as plt
import numpy as np
from previewAudio import previewAudio

SR = 40000

np.seterr(all='raise')

def simString(
    t_max = 2, pluck_height = .003, 
    pluck_pos = .5, 
    # pluck_pos = 1 / np.pi, 
    sr = SR, n_string_markers = 50, 
    spring_stiff = 921, string_density = .001, 
    string_stretch = 0.04, 
    wood_mass = 4, wood_damp = 0, 
    do_render = True, 
):
    time_step = 1 / sr
    marker_step = 1 / (n_string_markers + 1)
    rest_length = marker_step * (1 - string_stretch)
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
    if do_render:
        render(markers, n_markers, 0, pluck_height * 2)
    for clock in range(clock_max):
        t = clock * time_step
        markers[0,  0] = 0
        markers[-1, 0] = 1
        markers[0,  1] = wood_y
        markers[-1, 1] = wood_y
        displace = markers[1:] - markers[:-1]
        norm_displace = np.linalg.norm(displace, axis=1)
        unit_displace = np.divide(displace, norm_displace[:, np.newaxis])
        adj_displace = np.multiply(unit_displace, np.maximum(
            0, norm_displace - rest_length, 
        )[:, np.newaxis])
        wood_force = (adj_displace[0, 1] - adj_displace[-1, 1]) * (spring_stiff / marker_step)
        wood_velocity += wood_force * (time_step / wood_mass)
        wood_velocity += (- wood_damp * wood_velocity) * (time_step / wood_mass)
        wood_y += wood_velocity * time_step
        forces = (adj_displace[1:] - adj_displace[:-1]) * (spring_stiff / marker_step)
        markers_velocity += forces * (time_step / marker_mass)
        markers[1:-1] += markers_velocity * time_step
        audio[clock] = sampleAudio(wood_force)
        if do_render:
            render(markers, n_markers, t, pluck_height * 2)
    return audio

# time steps prt frame
TSPF = 1
# TSPF = 300
tspf_progress = 0
def render(markers, n_markers, t, y_ceil, force = False):
    global tspf_progress
    if not force:
        tspf_progress -= 1
        if tspf_progress <= 0:
            tspf_progress = TSPF
        else:
            return
    plt.clf()
    plt.plot(markers[:, 0], markers[:, 1])
    plt.plot(markers[:, 0], markers[:, 1], 'o')
    # plt.axis('equal')
    plt.axis([-.05, 1.05, -y_ceil, y_ceil])
    plt.title(f't = {t}')
    plt.draw()
    plt.pause(.15)

def sampleAudio(wood_force):
    return wood_force

audio = simString(do_render=1)
print('sim ok...')
previewAudio(audio, SR)
# plt.plot(audio)
# plt.show()

'''
4 kg:
47 cm -> 49 cm
921 N! 
0.41g
desity: .0002 kg / m
'''
