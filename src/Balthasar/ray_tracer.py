import re
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from manim import *
from obspy import read_inventory
from obspy.geodetics import gps2dist_azimuth
import ray_tracer_cpp as rtc

# === НАСТРОЙКИ ФАЙЛОВ ===
FCTOPT_PATH = "/home/ksukskss/Projects/NERV/data/FCTs_20100227_063411_NEAR_COAST_OF_CENTRAL_CHILE/fctoptsource_20100227_063411_NEAR_COAST_OF_CENTRAL_CHILE"
SEED_PATH = "/home/ksukskss/Projects/NERV/data/chil10058.seed"
manim_stations_dict = {}


#парсинг землетрясения из fctopt
def extract_earthquake_coords(fctopt_file):    
    with open(fctopt_file, "r") as f:
        first_line = f.readline().strip()    
    match = re.search(
        r"([-+]?\d*\.\d+|\b[-+]?\d+\b)\s+([-+]?\d*\.\d+|\b[-+]?\d+\b)$",
        first_line,
    )
    if not match:
        raise ValueError(
            f"Не удалось распарсить координаты из строки: {first_line}"
        )

    lat = float(match.group(1))
    lon = float(match.group(2))
    return lat, lon

#парсинг SEED + расчет эпицент растояний и азимутов
def get_stations_data(seed_file, eq_lat, eq_lon):
    inv = read_inventory(seed_file)

    stations_list = []
    seen_stations = set()

    for network in inv:
        for station in network:
            st_code = station.code
            if st_code in seen_stations:
                continue
            seen_stations.add(st_code)

            st_lat = station.latitude
            st_lon = station.longitude
            dist_m, azimuth, back_azimuth = gps2dist_azimuth(
                eq_lat, eq_lon, st_lat, st_lon
            )
            dist_deg = np.degrees(dist_m / 6371000.0)
            signed_azimuth = azimuth if azimuth <= 180 else azimuth - 360

            stations_list.append(
                {
                    "code": st_code,
                    "dist_deg": round(dist_deg, 1),
                    "azimuth": round(signed_azimuth, 1),
                }
            )

    return stations_list

#офрмление словаря для манима
def filter_and_format_stations(stations_raw):
    stations_raw.sort(key=lambda x: x["azimuth"])
    bins = np.linspace(-180, 180, 13)
    selected_stations = []

    for i in range(len(bins) - 1):
        low, high = bins[i], bins[i + 1]
        in_bin = [s for s in stations_raw if low <= s["azimuth"] < high]        
        # теневая зона 
        valid_in_bin = [s for s in in_bin if not (100.0 <= s["dist_deg"] <= 150.0)]
        if valid_in_bin:
            selected_stations.append(valid_in_bin[0])

    for st in selected_stations:
        label = f"{st['code']}"
        sign = 1 if st["azimuth"] >= 0 else -1
        manim_stations_dict[label] = (st["dist_deg"] * sign, "GREEN")

    #добавка фиктивных станция для красоты
    has_distant_station = any(150 <= abs(data[0]) <= 180 for data in manim_stations_dict.values())
    if not has_distant_station:
        manim_stations_dict["FICT_1"] = (150.0, "GREY")
        manim_stations_dict["FICT_2"] = (160.0, "GREY")
        manim_stations_dict["FICT_3"] = (-175.0, "GREY")
    return manim_stations_dict

#лог станций
eq_lat, eq_lon = extract_earthquake_coords(FCTOPT_PATH)
raw_data = get_stations_data(SEED_PATH, eq_lat, eq_lon)
final_stations = filter_and_format_stations(raw_data)
print("        stations = {")
for key, val in final_stations.items():

    print(f'            "{key}": ({val[0]}, {val[1]}),')
print("        }")



R_EARTH = 6371.0

# разбиение модели на слои 
def build_layers(csv_path):
    df = pd.read_csv(csv_path, header=None, names=["depth", "radius", "vp", "vs"])
    layers_p = []
    layers_s = []
    for i in range(len(df)-1):
        r1, vp1 = float(df.iloc[i]['radius']), float(df.iloc[i]['vp'])
        r2, vp2 = float(df.iloc[i+1]['radius']), float(df.iloc[i+1]['vp'])
        
        vs1 = float(df.iloc[i]['vs'])
        vs2 = float(df.iloc[i+1]['vs'])
        
        if r2 == 0.0: r2 = 1e-4
        if r1 == r2: continue
        
        #передача в С++ 
        layers_p.append(rtc.Layer(r1, vp1, r2, vp2))
        layers_s.append(rtc.Layer(r1, vs1, r2, vs2))

    return layers_p, layers_s

#ВИЗУЛ МАНИМ
class SeismicRayTracing(Scene):
    
    def construct(self):
        self.camera.background_color = WHITE
        SCALE = 3.0 / R_EARTH
        
        def pol2cart(r, delta_deg):
            # ВАЙБКОД ОГРАНИЧЕНИЕ ХЗ ДЛЯ ЧЕГО
            r = min(r, R_EARTH)
            delta_rad = np.radians(delta_deg)
            x = r * SCALE * np.sin(delta_rad)
            y = r * SCALE * np.cos(delta_rad)
            return np.array([x, y, 0])
        
        earth = Circle(radius=R_EARTH * SCALE, color=ORANGE, stroke_width=2)
        earth_fill = Circle(radius=R_EARTH * SCALE, color=ORANGE, fill_opacity=0.1, stroke_width=0)
        
        outer_core = Circle(radius=(R_EARTH - 2889) * SCALE, color=ORANGE, fill_opacity=0.2, stroke_width=2)
        DashedVMobject(outer_core, num_dashes=60)
        inner_core = Circle(radius=(R_EARTH - 5153) * SCALE, color=RED, fill_opacity=0.4, stroke_width=1)
        
        ev_depth = 150.0
        
        source_star = Star(color=RED, fill_opacity=1).scale(0.2).move_to(pol2cart(R_EARTH - ev_depth, 0))

        layers_p, layers_s = build_layers("/home/ksukskss/Projects/NERV/data/IASP91.csv")

        #набор для проверки
        stations = {
            "Япония (+45°)": (45.0, GREEN),
            "Япония (+60°)": (60.0, GREEN), "Япония (+90°)": (90.0, GREEN),
            "Япония (+100°)": (100.0, GREEN), "Япония (+155°)": (155.0, GREEN),
            "Япония (+170°)": (170.0, GREEN), "Япония (+180°)": (180.0, GREEN),
            "Япония (-20°)": (-20.0, GREEN), "Япония (-50°)": (-50.0, GREEN),
            "Япония (-67°)": (-67.0, GREEN), "Япония (-85°)": (-85.0, GREEN),
            "Япония (-150°)": (-150.0, GREEN),
            "Япония (-177°)": (-177.0, GREEN)
        }

        p_rays = []
        s_rays = []
        station_mobjects = []

        for st_name, (raw_dist, color) in manim_stations_dict.items():
            target_angle = raw_dist % 360
            is_right_side = target_angle <= 180
            shortest_dist = min(target_angle, 360 - target_angle)
            st_pos = pol2cart(R_EARTH, target_angle)
            st_marker = Triangle(color=GRAY, fill_opacity=1).scale(0.1).move_to(st_pos)
          
            st_marker.rotate(-np.radians(target_angle))
            direction_vector = normalize(st_marker.get_center())
            st_label = Text(st_name, font_size=16, color=GRAY)
            st_label.rotate(-np.radians(target_angle))
            st_label.rotate(PI / 2)
            if 180 < target_angle < 360:
                st_label.rotate(PI)
            st_label.next_to(st_marker, direction_vector, buff=0.15)
            station_group = VGroup(st_marker, st_label)
            station_mobjects.append(station_group)
            ray_results_p = rtc.find_all_takeoff_angles(shortest_dist, ev_depth, layers_p)
            for res in ray_results_p:
                traveled_dist = res.traveled_dist
                path_r = res.path.r           # Массив радиусов из C++
                path_delta = res.path.delta   # Массив углов из C++
                if abs(path_r[-1] - R_EARTH) > 5.0:
                    continue

                is_short_ray = traveled_dist <= 180
                multiplier = 1 if (is_right_side == is_short_ray) else -1
                points = [pol2cart(r, np.degrees(d) * multiplier) for r, d in zip(path_r, path_delta)]
                
                ray_line = VMobject().set_points_as_corners(points).set_color(YELLOW).set_stroke(width=2.5)
                p_rays.append(ray_line)

            # расчет лучей 
            ray_results_s = rtc.find_all_takeoff_angles(shortest_dist, ev_depth, layers_s)
            for res in ray_results_s:
                traveled_dist = res.traveled_dist
                path_r = res.path.r
                path_delta = res.path.delta

                if abs(path_r[-1] - R_EARTH) > 5.0:
                    continue
                is_short_ray = traveled_dist <= 180
                multiplier = 1 if (is_right_side == is_short_ray) else -1
                points = [pol2cart(r, np.degrees(d) * multiplier) for r, d in zip(path_r, path_delta)]
                solid_ray_s = VMobject().set_points_as_corners(points)
                ray_line_s = DashedVMobject(solid_ray_s, num_dashes=25).set_color(BLACK).set_stroke(width=3.0)
                s_rays.append(ray_line_s)

        source_star.set_z_index(10)
        for st in station_mobjects:
            st.set_z_index(10)
        self.add(earth, earth_fill, outer_core, inner_core, source_star, *station_mobjects)
        self.wait(0.5)

        p_animations = [Create(ray, rate_func=linear).set_run_time(4.0) for ray in p_rays]
        s_animations = [Create(ray, rate_func=linear).set_run_time(4.0) for ray in s_rays]
        
        # БЕЗОПАСНЫЙ ЗАПУСК
        if p_animations or s_animations:
            self.play(
                *p_animations,
                *s_animations
            )
        else:
            print("ВНИМАНИЕ: Ни одного луча не было сгенерировано!")
            
        self.wait(3)
