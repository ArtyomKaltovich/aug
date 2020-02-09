import itertools
import random
import matplotlib.pyplot as plt

from aug.qc.gc import gc_plot
from aug.seq.seq import gc_rate


def gen_data():
    bunch_size = 100
    gc_percent_min = 40
    gc_percent_max = 60
    # gen uniform data
    data1 = []
    for _, gc_percent in zip(range(bunch_size), itertools.cycle(range(gc_percent_min, gc_percent_max))):
        at_percent = 100 - gc_percent
        read = "".join(random.choices("CG", k=gc_percent)) + "".join(random.choices("AT", k=at_percent))
        print(gc_rate(read), read)
        data1.append(read)
    #data2 = ["".join(random.choices("ACGT", k=100)) for _ in range(bunch_size)]
    data2 = data1
    return data1, data2


if __name__ == '__main__':
    data1, data2 = gen_data()
    fig, (ax1, ax2) = plt.subplots(1, 2)
    gc_plot(data1, ax=ax1)
    gc_plot(data2, ax=ax2)
    plt.show()
