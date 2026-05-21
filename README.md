NERV — **N**umerical **E**arth **R**ay **V**isualizer. 

Sructure
nerv-seismo/
│
├── data/
│   ├── fctoptsource_20100227_063411_NEAR_COAST_OF_CENTRAL_CHILE   
│   └── iasp91.csv                                                  
│
├── src/
│   ├── __init__.py
│   │
│   ├── magi/
│   │   ├── __init__.py
│   │   │
│   │   ├── caspar/          # Задача 1: Географический домен
│   │   │   ├── __init__.py
│   │   │   └── geo_calc.py    # Сферическая геометрия, расстояния, азимуты
│   │   │
│   │   ├── balthasar/         # Задача 2: Главный фокус №1 (Shooting Ray-Tracer)
│   │   │   ├── __init__.py
│   │   │   ├── earth_layers.py# Дискретизация IASP91 на градиентные слои
│   │   │   └── shooter.py     # Алгоритм пристрелки и вычисление интегралов пути
│   │   │
│   │   └── melchior/            # Задача 3: Главный фокус №2 (Convolution & Synthesis)
│   │       ├── __init__.py
│   │       ├── radiation.py   # Вспомогательный векторный расчет механизма очага
│   │       ├── rotation.py    # Вспомогательные матрицы поворота волн на приемнике
│   │       └── convolve.py    # Мастер-цикл дискретной свертки во временной области
│   │
│   ├── utils/
│   │   ├── __init__.py
│   │   └── scardec_parser.py  # Парсер заголовков и STF-массива из файла
│   │
│   └── main.py                 # Запуск терминала NERV (Центральный оркестратор)
│
│
└── README.md

For beginning, need to compare 3 major elements of data. (I think you have been created the project directory =))
First needs find and download file with earthquakes parameters and STF. 
For Chile 2010 event:

'''
mkdir data/
cd data/
curl -LOC - http://scardec.projects.sismo.ipgp.fr/arch/sourcefunction_archive_152924350.tar.gz && tar -xvzf sourcefunction_archive_152924350.tar.gz 
'''

Next download the SEED or miniSEED file.
GEOSCOPE network for Chile:

'''
curl -LOC - http://geoscope.ipgp.fr/seismes/SEED/G/2010/chil10058.seed.gz && gunzip chil10058.seed.gz
'''

