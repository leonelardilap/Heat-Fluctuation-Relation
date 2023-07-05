# Heat exchange fluctuation relation 🔥

<p align="center">
  <img align="left" width="400" height="400"src="simulation/animation/ensembles.gif" alt="animated" />
</p>

Here we introduce the simulation of the transition from a microcanonical to a canonical ensembles in the classical regime. the simulation it's done through molecular dynamics simulations for the microcanonical ensemble of a single harmonic oscillator evolving into a canonical ensemble in presence of a thermal bath.

## Compile and run the simulation

Compile the simulation code with the first line and then with the second line execute it

```bash
g++ -o simulation micro2cano.cpp
./simulation
```

## Output files

- mean_energy.dat
- microcanonical_states.dat
- canonical_states.dat
- transfer_heat.dat

## Python enviroment

To create the python environment used the .yml file. Write into the console

```bash
conda env create -f environment.yml
```

[//]: # "To export the created conda environment into the .yml file, with your conda environment activated, run the following command to generate dependency yaml file"
[//]: # "conda env export > environment.yml"
