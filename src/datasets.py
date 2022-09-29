import pathlib
import re
import pandas as pd
import numpy as np
from eaps import KeplerianElements

# ----------------------------------------------------------------------------------------
# Path Utilities
# ----------------------------------------------------------------------------------------
def data_directory(*args) -> pathlib.Path:
    pth = pathlib.Path("__file__").parent.resolve().parent.resolve() / "data"

    for arg in args:
        pth /= arg

    return pth


# ----------------------------------------------------------------------------------------
# Data set
# ----------------------------------------------------------------------------------------
class Dataset:
    def count(self) -> int:
        return self.data.shape[0]

    @classmethod
    def load_all(cls, dir: pathlib.Path) -> list:
        datasets = []

        if not dir.is_dir():
            return datasets

        for file in dir.iterdir():
            if re.match(cls.REGEX_FILTER, file.name):
                datasets.append(cls(file.resolve()))

        return datasets


class StateDataset(Dataset):
    COLUMN_NAMES = ["t", "x", "y", "z", "vx", "vy", "vz"]
    REGEX_FILTER = re.compile(r"id\d{6}-XV.csv")

    def __init__(self, path: pathlib.Path):
        self.path = path
        self.name = path.name
        self.data = pd.read_csv(self.path, skiprows=1, names=StateDataset.COLUMN_NAMES)

    def times_and_states(self) -> np.ndarray:
        out = np.zeros((6, self.count()))
        out[0, :] = self.data.x
        out[1, :] = self.data.y
        out[2, :] = self.data.z
        out[3, :] = self.data.vx
        out[4, :] = self.data.vy
        out[5, :] = self.data.vz
        return (self.data.t.to_numpy(), out)

    def output_path(self) -> pathlib.Path:
        parent = self.path.parent
        name = self.path.stem.strip("-XV") + "-EL" + self.path.suffix
        return parent / name


class ElementDataset(Dataset):
    COLUMN_NAMES = [
        "t",
        "a",
        "e",
        "inclination",
        "lon_asc_node",
        "arg_peri",
        "true_anom",
    ]
    REGEX_FILTER = re.compile(r"id\d{6}-EL.csv")

    def __init__(self, path: pathlib.Path):
        self.path = path
        self.name = path.name
        self.data = pd.read_csv(self.path, skiprows=1, names=ElementDataset.COLUMN_NAMES)

    def times_and_elements(self) -> np.ndarray:
        out = KeplerianElements(
            self.data.a,
            ecc=self.data.e,
            inc=self.data.inclination,
            raan=self.data.lon_asc_node,
            aop=self.data.arg_peri,
            ta=self.data.true_anom
        )
        return (self.data.t.to_numpy(), out)

    def output_path(self) -> pathlib.Path:
        parent = self.path.parent.resolve().parent.resolve() / "plots"
        name = self.path.stem.strip("-XV") + "-EL" + ".png"
        return parent / name