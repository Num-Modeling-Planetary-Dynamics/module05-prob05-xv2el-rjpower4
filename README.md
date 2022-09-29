[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-c66648af7eb3fe8bc4f294546bfd86ef473780cde1dea487d3c4ff354943c9ae.svg)](https://classroom.github.com/online_ide?assignment_repo_id=8589798&assignment_repo_type=AssignmentRepo)
| EAPS 591 - Numerical Modeling of Planetary Orbits | Fall 2022 | [Rolfe Power] |
| ----------------------------- | --------- | ------------------ |

# Usage

You'll need python with the following packages

- numpy
- matplotlib
- pandas

To create the output csv files, simply `cd` into the `src` directory and run

```
$ python ./conversion.py
```

To create the plots, make sure you've run the `conversion.py` script and then while in the `src` directory, run

```
$ python ./make_plots.py
```