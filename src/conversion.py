import pandas as pd
import eaps
from datasets import StateDataset, ElementDataset, data_directory

# ----------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------
GM_SUN = 2.95912208285590931905e-04 # AU^3 / DAY^2

# ----------------------------------------------------------------------------------------
# Conversion
# ----------------------------------------------------------------------------------------
def process_state_dataset(ds: StateDataset):
    (t, ss) = ds.times_and_states()
    elements = eaps.KeplerianElements.from_states(1.0, ss)
    return pd.DataFrame(
        data={
            "t": t,
            "a": elements.semi_major_axis(),
            "e": elements.eccentricity(),
            "inclination": elements.inclination(),
            "lon_asc_node": elements.right_ascension(),
            "arg_peri": elements.argument_of_periapsis(),
            "true_anom": elements.true_anomaly(),
        }
    )

def process_element_dataset(ds: ElementDataset):
    (t, els) = ds.times_and_elements()
    states = els.to_inertial_states(GM_SUN)

    return pd.DataFrame(
        data = {
             "t": t,
             "xh": states[0],
             "yh": states[1],
             "zh": states[2],
             "vxh": states[3],
             "vyh": states[4],
             "vz": states[5]
        }
    )


# ----------------------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------------------
def main():
    datasets = ElementDataset.load_all(data_directory())
    results = [process_element_dataset(ds) for ds in datasets]
    for (ds, res) in zip(datasets, results):
        output_path = ds.data_output_path()
        res.to_csv(output_path)

if __name__ == "__main__":
    main()