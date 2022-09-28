import pandas as pd
import eaps
from datasets import StateDataset, data_directory

# ----------------------------------------------------------------------------------------
# Conversion
# ----------------------------------------------------------------------------------------
def process_dataset(ds: StateDataset):
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


# ----------------------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------------------
def main():
    datasets = StateDataset.load_all(data_directory())
    results = [process_dataset(ds) for ds in datasets]

    for (ds, res) in zip(datasets, results):
        res.to_csv(ds.output_path())

if __name__ == "__main__":
    main()