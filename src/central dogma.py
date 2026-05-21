# src/central_dogma.py
import os
from Balthasar.geometry import (
    parse_scardec_source, extract_station_metadata, calculate_epicentral_angles
)
from Balthasar.ray_tracer import BalthasarRayTracer

def run_pipeline():
    # Настраиваем точные пути к вашим файлам
    SCARDEC_FILE = "data/FCTs_20100227_063411_NEAR_COAST_OF_CENTRAL_CHILE/fctoptsource_20100227_063411_NEAR_COAST_OF_CENTRAL_CHILE"
    SEED_FILE = "data/chil10058.seed"  # Ваш SEED файл в корне или data/raw/
    IASP91_CSV = "data/IASP91.csv"      # Ваша модель скоростей

    print("[NERV] Запуск геодезического анализа...")
    
    # 1. Парсим источник
    source_meta = parse_scardec_source(SCARDEC_FILE)
    
    # 2. Парсим метаданные станции SCZ из вашего SEED-файла
    station_meta = extract_station_metadata(SEED_FILE, network="G", station="SCZ")
    
    # 3. Считаем углы большого круга
    geo_metrics = calculate_epicentral_angles(source_meta, station_meta)

    print(f"\n--- ГЕОМЕТРИЧЕСКИЙ КОНТЕКСТ БАЛЬТАЗАРА ---")
    print(f"Эпицентр (SCARDEC): Lat {source_meta['lat']}°, Lon {source_meta['lon']}°")
    print(f"Станция (SEED): SCZ | Lat {station_meta['lat']}°, Lon {station_meta['lon']}°")
    print(f"Глубина очага:       {source_meta['depth']} км")
    print(f"Истинное расстояние: {geo_metrics['delta_deg']:.4f}° ({geo_metrics['delta_rad']:.4f} рад)")

    # 4. Запуск трассировщика на вашей модели IASP91
    print("\n[Balthazar] Запуск численного трассировщика (Метод пристрелки)...")
    tracer = BalthasarRayTracer(IASP91_CSV)
    
    try:
        ray_results = tracer.find_ray_for_station(
            target_delta_rad=geo_metrics["delta_rad"],
            source_depth=source_meta["depth"]
        )
        
        print(f"\n--- СХОДИМОСТЬ ТРАССИРОВКИ ЛУЧЕЙ БАЛЬТАЗАРА ---")
        print(f"Тип волны:               Прямая P-волна")
        print(f"Расчетное время пробега: {ray_results['travel_time']:.2f} сек")
        print(f"Угол выхода (i_h):       {ray_results['takeoff_angle']:.2f}° (от вертикали вниз)")
        print(f"Угол прихода (i_rec):     {ray_results['incidence_angle']:.2f}° (к вертикали вверх)")
        print(f"Параметр луча (p):       {ray_results['ray_parameter']:.4f} с/км")
        
    except Exception as e:
        print(f"[Ошибка Бальтазара]: {e}")

if __name__ == "__main__":
    run_pipeline()
