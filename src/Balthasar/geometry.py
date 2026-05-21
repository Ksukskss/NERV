# src/magi_balthazar/geometry.py
import numpy as np
from obspy import read_inventory

def parse_scardec_source(scardec_path):
    """
    Извлекает метаданные очага из файла SCARDEC.
    """
    with open(scardec_path, 'r') as f:
        lines = f.readlines()
    
    tokens_l1 = lines[0].split()
    src_lat = float(tokens_l1[6])
    src_lon = float(tokens_l1[7])
    
    tokens_l2 = lines[1].split()
    src_depth = float(tokens_l2[0])
    
    return {"lat": src_lat, "lon": src_lon, "depth": src_depth}

def extract_station_metadata(seed_path, network, station):
    """
    Извлекает координаты станции напрямую из вашего бинарного SEED-файла.
    """
    # ObsPy автоматически парсит бинарный формат SEED и строит инвентарь
    inv = read_inventory(seed_path, format="SEED")
    sta_obj = inv.select(network=network, station=station)[0][0]
    
    return {
        "lat": sta_obj.latitude,
        "lon": sta_obj.longitude,
        "elevation": sta_obj.elevation / 1000.0  # Из метров в км
    }

def calculate_epicentral_angles(src_meta, sta_meta):
    """
    Вычисляет угловое расстояние (Delta) и азимуты.
    """
    from obspy.geodetics import gps2dist_azimuth, locations2degrees
    
    delta_deg = locations2degrees(src_meta["lat"], src_meta["lon"], 
                                   sta_meta["lat"], sta_meta["lon"])
    
    dist_m, azimuth, back_azimuth = gps2dist_azimuth(src_meta["lat"], src_meta["lon"], 
                                                     sta_meta["lat"], sta_meta["lon"])
    
    return {
        "delta_deg": delta_deg,
        "delta_rad": np.radians(delta_deg),
        "azimuth": azimuth,
        "back_azimuth": back_azimuth
    }