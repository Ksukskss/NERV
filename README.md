# NERV — Numerical Earthquake Ray Visualizer

### Ray-tracer is ready (Balthasar)!
Presintation: https://docs.google.com/presentation/d/1LTbex2gu9NkrWASXG-jDrI-VzhO9Dwuiu2QHmrWR9lk/edit?usp=sharing

A seismic ray-tracing and visualization tool that computes and animates the paths of P and S waves traveling through a 1-D Earth (IASP91 model). Given an earthquake hypocenter and a set of recording stations, it solves the two-point boundary-value problem — finding the exact take-off angle that connects the source to each station — and renders the resulting curved ray paths as a Manim animation. The numerically intensive ray tracing runs in C++ (exposed to Python via pybind11), while data parsing, geodesy, and visualization are handled in Python.

## Data Installation

To run this seismic ray-tracing script successfully, you need to prepare and place several external data files (earthquake parameters, station metadata, and 1D velocity models).
### 1. Chile 2010 Event Parameters and STF

```bash
mkdir -p data/
cd data/
curl -LOC - http://scardec.projects.sismo.ipgp.fr/arch/sourcefunction_archive_152924350.tar.gz && tar -xvzf sourcefunction_archive_152924350.tar.gz 
```

### 2. SEED / miniSEED Data

Next, download the SEED or miniSEED file. GEOSCOPE network for Chile:

```bash
curl -LOC - http://geoscope.ipgp.fr/seismes/SEED/G/2010/chil10058.seed.gz && gunzip chil10058.seed.gz
```

### 3. IASP91 Earth Model

Download the IASP91 model of Earth:

```bash
curl -LOC https://ds.iris.edu/spudservice/data/9991804 -o iasp91.csv
```

### Required Files Structure

Ensure your local data directory matches the paths referenced in the script or update the paths in the `# === FILE SETTINGS ===` section of the code.

```
Projects/
└── NERV/
    └── data/
        ├── FCTs_20100227_063411_NEAR_COAST_OF_CENTRAL_CHILE/
        │   └── fctopt_20100227_063411_NEAR_COAST_OF_CENTRAL_CHILE
        ├── chil10058.seed
        └── IASP91.csv
```

## Getting Started
Cloning the Repository:

```bash
git clone https://github.com/Ksukskss/NERV.git
cd NERV
```

Create a virtual environment:
```bash
uv venv
```

Install the required Python packages:
```bash
pip install numpy pandas scipy obspy manim
```
**Important:** The script relies on a custom C++ module called `ray_tracer_cpp`. Make sure you have compiled and installed this module into your Python environment

Running the Script:
```bash
uv run manim -pqh ray_tracer.py SeismicRayTracing
```
## Concept

The architectural solution and task division logic in the project are based on the mental model of the MAGI supercomputer (from the Neon Genesis Evangelion universe). The system has three independent computational modules that solve three different tasks:    
#### MAGI-1: BALTAZAR (seismic ray tracer)
>    Responsible for wave kinematics and dynamics. The module performs ray tracing and simulates the transfer of energy from the hypocentre to the Earth's surface through the heterogeneous mantle
>    
#### MAGI-2: CASPER (fault geometry and dynamics)  
#### MAGI-3: MELCHIOR 
>The module calculates (synthetic) seismograms using the kinematic parameters of the rays from Baltazar and the fault model from Casper.

## Structure

```text
NERV/
│
├── data/
│   ├── fctoptsource_20100227_063411_NEAR_COAST_OF_CENTRAL_CHILE
│   ├── iasp91.csv
│   └── chil10058.seed
│
├── src/
│   ├── __init__.py
│   │
│   ├── balthazar/          # THE PROPAGATION PATH ENGINE
│   │   ├── __init__.py
│   │   ├── ray_tracer.cpp       #ray tracer on c++
│   │   └── ray_tracer.py        #python 
│   │
│   ├── casper/             # THE COORDINATE & ROTATION ENGINE
│   │   ├── __init__.py
│   │   ├── source_radiation.py  # Fault vectors and radiation patterns (P, SV, SH)
│   │   └── receiver_rotation.py # ZRT to ZNE projection matrices
│   │
│   ├── melchior/           # THE SIGNAL PROCESSING KERNEL
│   │   ├── __init__.py
│   │   └── convolution.py       # Direct time-domain sliding convolution loop
│   │
│   └── central_dogma.py         # Main orchestrator running the pipeline loops
│
└── README.md
```


