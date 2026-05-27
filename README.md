# NERV — Numerical Earth Ray Visualizer

https://docs.google.com/presentation/d/1LTbex2gu9NkrWASXG-jDrI-VzhO9Dwuiu2QHmrWR9lk/edit?usp=sharing

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
│   │   ├── geometry.py          # Spherical distance & azimuth calculations
│   │   └── ray_tracer.py        # Shooting/bisection 1D layer tracer
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

## Data Setup

To begin, you need to compare 3 major elements of data. First, find and download the file with earthquake parameters and STF.

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
