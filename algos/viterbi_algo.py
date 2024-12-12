import numpy as np


def viterbi(obs, trans, emission, initialProb, emission_key):
    """
    :param obs: np array, in sequence order
    :param trans: double numpy array
    :param emission: double numpy array, each row represents a state
    :param initialProb: arr, 0 for header(single), 1 for time delta, 2 for on, 3 for note, 4 for velocity,
    5 for intergenic(single)
    :param emission_key:
    :return:
    """
    t = len(obs)
    k = trans.shape[0]

    viterbi_table = np.empty((k, t))

    # The first column
    viterbi_table[:, 0] = np.log(initialProb.T * emission[:, emission_key.index(obs[0])])

    for i in range(1, t):
        viterbi_table[0, i] = viterbi_table[0, i - 1] + np.log(trans[0, 0]) + np.log(
            emission[0, emission_key.index(obs[i])])
        viterbi_table[1, i] = max(
            viterbi_table[0, i - 1] + np.log(trans[0, 1]) + np.log(emission[1, emission_key.index(obs[i])]),
            viterbi_table[5, i - 1] + np.log(trans[5, 1]) + np.log(emission[1, emission_key.index(obs[i])])
        )
        viterbi_table[2, i] = viterbi_table[1, i - 4] + np.log(
            trans[1, 2] + np.log(emission[2, emission_key.index(obs[i])]))
        viterbi_table[3, i] = viterbi_table[2, i - 4] + np.log(
            trans[2, 3] + np.log(emission[3, emission_key.index(obs[i])]))
        viterbi_table[4, i] = viterbi_table[3, i - 4] + np.log(
            trans[3, 4] + np.log(emission[4, emission_key.index(obs[i])]))
        viterbi_table[5, i] = max(
            viterbi_table[5, i - 1] + np.log(trans[5, 5] + np.log(emission[5, emission_key.index(obs[i])])),
            viterbi_table[4, i - 4] + np.log(trans[4, 5] + np.log(emission[5, emission_key.index(obs[i])])),
        )

    # backtracking starts
    path = np.zeros(t)
    path[t - 1] = np.argmax(viterbi_table[:, -1])

    i = t - 1
    while i > 0:
        if path[i] == 0:
            path[i - 1] = 0
            i -= 1
        elif path[i] == 1:
            path[i - 1] = 1
            path[i - 2] = 1
            path[i - 3] = 1
            if viterbi_table[0, i - 1] + np.log(trans[0, 1]) + np.log(emission[1, emission_key.index(obs[i])]) > \
                    viterbi_table[5, i - 1] + np.log(trans[5, 1]) + np.log(emission[1, emission_key.index(obs[i])]):
                path[i - 4] = 0
            else:
                path[i - 4] = 0
            i -= 4

        elif path[i] == 2:
            path[i - 1] = 2
            path[i - 2] = 2
            path[i - 3] = 2
            path[i - 4] = 1
            i -= 4
        elif path[i] == 3:
            path[i - 1] = 3
            path[i - 2] = 3
            path[i - 3] = 3
            path[i - 4] = 2
            i -= 4
        elif path[i] == 4:
            path[i - 1] = 4
            path[i - 2] = 4
            path[i - 3] = 4
            path[i - 4] = 3
            i -= 4
        elif path[i] == 5:
            if viterbi_table[5, i - 1] + np.log(trans[5, 5] + np.log(emission[5, emission_key.index(obs[i])]))> \
                    viterbi_table[4, i - 4] + np.log(trans[4, 5] + np.log(emission[5, emission_key.index(obs[i])])):
                path[i-1] = 5
            else:
                path[i-1] = 4
            i -= 1
    return path, viterbi_table

def main():
    
