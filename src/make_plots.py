from datasets import RoundTripStateDataset, StateDataset, data_directory
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import re

def plot_dataset(title, ds: RoundTripStateDataset, sds: StateDataset):
    
    (_, s1) = ds.times_and_states()
    (_, s2) = sds.times_and_states()

    dr = np.abs(norm(s1[:3], axis=0) - norm(s2[:3], axis=0))
    dv = np.abs(norm(s1[3:], axis=0) - norm(s2[3:], axis=0))

    fig, axs = plt.subplots(2, sharex=True)
    fig.set_size_inches(10, 8)
    fig.suptitle(title, fontsize=20)
    axs[0].plot(ds.data.t, dr)
    axs[0].set_ylabel("Position Magnitude Error", fontsize=18)
    axs[0].tick_params('both', labelsize=18)
    axs[0].set_yscale("log")

    axs[1].plot(ds.data.t, dv)
    axs[1].set_ylabel("Velocity Magitude Error", fontsize=18)
    axs[1].set_xlabel("Time", fontsize=18)
    axs[1].tick_params('both', labelsize=18)
    axs[1].set_yscale("log")

    return (fig, axs)


def get_pairs() -> dict:
    dd = data_directory()
    reg = re.compile(r"(id\d{6})-XV.csv")
    pairs = {}
    for k in dd.iterdir():
        m = re.match(reg, k.name)
        if m:
            root = m.groups()[0]
            pairs[root] = (
                (dd / (root + "-XV.csv")),
                (dd / (root + "-EL-TO-XV.csv"))
            )
    return pairs

# ----------------------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------------------
def main():
    pairs = get_pairs()
    for key in pairs.keys():
        rtsds = RoundTripStateDataset(pairs[key][1])
        sds = StateDataset(pairs[key][0])
        f, a = plot_dataset(key, rtsds, sds)
        f.savefig(rtsds.output_path(), dpi=300)

if __name__ == "__main__":
    main()
# if __name__ == "__main__":
#     datasets = main()
#     for ds in datasets:
#         f, a = plot_dataset(ds)
#         f.savefig(ds.output_path(), dpi=300)