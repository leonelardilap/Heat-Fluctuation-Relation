# Heat exchange fluctuation relation 🔥

<p>
  <img align="right" width="400" height="400"src="simulation/animation/ensembles.gif" alt="animated" />
</p>

In this repo through molecular dynamics techniques, in the classical regime, is simulated the microcanonical ensemble of a single harmonic oscillator evolving to a canonical one in presence of a thermal. with the results of simulations is validated numerically a fluctuation relation for the heat exchange between a microcanonical and a canonical ensembles.

## Compile and run the simulation

Compile the code in the simulation folder with the first line and then with the second line execute it

```bash
g++ -o heatup ensemble.cpp
./heatup
```

## Python enviroment

To create the python environment used the .yml file. Write into the console

```bash
conda env create -f environment.yml
```

[//]: # "To export the created conda environment into the .yml file, with your conda environment activated, run the following command to generate dependency yaml file"
[//]: # "conda env export > environment.yml"
