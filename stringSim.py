'''
'''
from matplotlib import pyplot as plt
import numpy as np
import scipy.signal
from previewAudio import previewAudio
from jdt import Jdt
from yin import yin

SR = 100000

np.seterr(all='raise')

def simString(
    t_max = 2, pluck_height = .005, 
    # pluck_pos = .5, 
    pluck_pos = 1 / np.pi, 
    sr = SR, n_string_markers = 90, 
    spring_stiff = 921, string_density = .001, 
    string_stretch = 0.18, 
    wood_mass = 4, wood_damp = 0, 
    unbendability = 1, 
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
    wood_y = np.zeros((2, ))
    wood_velocity = np.zeros((2, ))
    clock_max = t_max * sr
    audio = np.zeros((clock_max, ))
    if do_render:
        render(markers, n_markers, 0, pluck_height * 2)
    with Jdt(clock_max, UPP = 2048) as j:
        for clock in range(clock_max):
            t = clock * time_step
            markers[0,  0] = 0
            markers[-1, 0] = 1
            markers[0,  1] = wood_y[0]
            markers[-1, 1] = wood_y[1]
            displace = markers[1:] - markers[:-1]
            norm_displace = np.linalg.norm(displace, axis=1)
            unit_displace = np.divide(displace, norm_displace[:, np.newaxis])
            adj_displace = np.multiply(unit_displace, np.maximum(
                0, norm_displace - rest_length, 
            )[:, np.newaxis])
            wood_velocity += spring_stiff / marker_step * np.array([
                adj_displace[0, 1], - adj_displace[-1, 1]
            ]) * (time_step / wood_mass)
            # wood_velocity += (- wood_damp * wood_velocity) * (time_step / wood_mass)
            wood_velocity *= np.exp(- wood_damp * time_step)
            wood_y += wood_velocity * time_step
            forces = (adj_displace[1:] - adj_displace[:-1]) * (
                spring_stiff / marker_step
            ) + (unit_displace[1:] - unit_displace[:-1]) * unbendability
            markers_velocity += forces * (time_step / marker_mass)
            markers[1:-1] += markers_velocity * time_step
            # audio[clock] = wood_force
            audio[clock] = wood_y[0]
            # audio[clock] = adj_displace[0, 1]
            if do_render:
                render(markers, n_markers, t, pluck_height * 2)
            j.acc()
    return audio

# time steps prt frame
# TSPF = 1
TSPF = 550
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
    plt.title(f't = {format(t, ".03f")}')
    plt.draw()
    plt.pause(.15)

def sampleAudio(wood_force):
    return wood_force

def contrast(audio, x = 404, ax = plt):
    hann = scipy.signal.get_window('hann', audio.size, True)
    spec = np.abs(np.fft.rfft(hann * audio))
    try:
        log_spec = np.log(spec)
    except FloatingPointError:
        np.seterr(all='warn')
        log_spec = np.log(spec)
        np.seterr(all='raise')
    ax.plot(log_spec)
    for i in range(45):
        plt.axvline(x * i, c='r')
    ax.axis([-100, 19000, -22, -9])

def studyUnbendability():
    trials = [0, 1, 5, 50, 200]
    fig, axes = plt.subplots(len(trials))
    for unbendability, ax in zip(trials, axes):
        print('unbendability:', unbendability)
        audio = simString(do_render=0, unbendability=unbendability)
        contrast(audio, ax=ax)
    fig.show()

studyUnbendability()

'''
4 kg:
47 cm -> 49 cm
921 N! 
0.41g
desity: .0002 kg / m
'''
