import re
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from manim import *
from obspy import read_inventory
from obspy.geodetics import gps2dist_azimuth

# === НАСТРОЙКИ ФАЙЛОВ ===
FCTOPT_PATH = "/home/ksukskss/Projects/NERV/data/FCTs_20100227_063411_NEAR_COAST_OF_CENTRAL_CHILE/fctoptsource_20100227_063411_NEAR_COAST_OF_CENTRAL_CHILE"
SEED_PATH = "/home/ksukskss/Projects/NERV/data/chil10058.seed"
manim_stations_dict = {}

def extract_earthquake_coords(fctopt_file):
    """Читает первую строчку fctopt и извлекает шир. и долг. землетрясения."""
    with open(fctopt_file, "r") as f:
        first_line = f.readline().strip()

    # Ищем вещественные числа в конце строки времени (например, -36.120 -72.900)
    # Формат первой строки обычно: YYYY MM DD HH MM SS.S LAT LON
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


def get_stations_data(seed_file, eq_lat, eq_lon):
    """Читает SEED инвентарь, вычисляет угловые расстояния (в градусах)

    и азимуты от очага до каждой уникальной станции.
    """
    # Загружаем метаданные станций из SEED-файла
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

            # gps2dist_azimuth возвращает: (дистанция_в_метрах, прямой_азимут, обратный_азимут)
            # Прямой азимут (Great Circle) — это угол от севера из эпицентра в сторону станции.
            dist_m, azimuth, back_azimuth = gps2dist_azimuth(
                eq_lat, eq_lon, st_lat, st_lon
            )

            # Переводим метры в градусы дуги большой окружности земли (1 градус ~ 111.195 км)
            # Или точнее, используя средний радиус Земли:
            dist_deg = np.degrees(dist_m / 6371000.0)

            # Для 1D проекции в Manim нам часто нужно знать "направление" луча.
            # Если мы хотим развернуть станции на плоскости от -180° до +180°
            # на основе азимута вылета волны:
            # Превращаем азимут [0, 360] в диапазон [-180, 180] для сортировки
            signed_azimuth = azimuth if azimuth <= 180 else azimuth - 360

            stations_list.append(
                {
                    "code": st_code,
                    "dist_deg": round(dist_deg, 1),
                    "azimuth": round(signed_azimuth, 1),
                }
            )

    return stations_list
def filter_and_format_stations(stations_raw):
    """Фильтрует станции, чтобы получить равномерное покрытие от -180 до 180 градусов,

    выбирая около 10-12 станций с разными расстояниями.
    """
    # Сортируем станции по значению знакового азимута (от -180 до +180)
    stations_raw.sort(key=lambda x: x["azimuth"])

    # Будем группировать станции по "корзинам" азимутов, чтобы гарантировать круговой охват
    # Разделим окружность на 12 секторов по 30 градусов
    bins = np.linspace(-180, 180, 25)
    selected_stations = []

    for i in range(len(bins) - 1):
        low, high = bins[i], bins[i + 1]
        # Ищем все станции, попавшие в этот азимутальный сектор
        in_bin = [s for s in stations_raw if low <= s["azimuth"] < high]
        valid_in_bin = [
            s
            for s in in_bin
            if not (100.0 <= abs(s["azimuth"]) <= 150.0)
        ]

        # Если в секторе остались подходящие станции вне теневой зоны, берем первую
        if valid_in_bin:
            selected_stations.append(valid_in_bin[0])

    # Если станций вышло чуть больше/меньше, можно взять срез или добавить условия,
    # но разбиение по корзинам автоматически даст равномерный красивый круг.

    # Форматируем в итоговый словарь строк для Manim
    
    for st in selected_stations:
        # Формируем имя вида "CODE (+45.2°)" или "CODE (-120.5°)"
        sign = "+" if st["azimuth"] >= 0 else ""
        label = f"{st['code']}"

        # Добавляем в словарь. dist_deg — эпицентральное расстояние, используемое Manim для угла
        manim_stations_dict[label] = (st["azimuth"], "GREEN")

    return manim_stations_dict
eq_lat, eq_lon = extract_earthquake_coords(FCTOPT_PATH)
# 2. Читаем SEED и считаем геодезию
raw_data = get_stations_data(SEED_PATH, eq_lat, eq_lon)
# 3. Фильтруем для кругового 1D охвата
final_stations = filter_and_format_stations(raw_data)
print("        stations = {")
for key, val in final_stations.items():

    print(f'            "{key}": ({val[0]}, {val[1]}),')
print("        }")
# ==========================================
# 1. МАТЕМАТИЧЕСКАЯ ЧАСТЬ (Аналитика)
# ==========================================
R_EARTH = 6371.0

def build_layers(csv_path):
    df = pd.read_csv(csv_path, header=None, names=["depth", "radius", "vp", "vs"])
    layers_p = []
    layers_s = []
    for i in range(len(df)-1):
        r1, vp1 = float(df.iloc[i]['radius']), float(df.iloc[i]['vp'])
        r2, vp2 = float(df.iloc[i+1]['radius']), float(df.iloc[i+1]['vp'])
        
        # Считываем также скорости S-волн
        vs1 = float(df.iloc[i]['vs'])
        vs2 = float(df.iloc[i+1]['vs'])
        
        if r2 == 0.0: r2 = 1e-4  # Защита от деления на ноль в центре Земли
        if r1 == r2: continue
        
        layers_p.append((r1, vp1, r2, vp2))
        layers_s.append((r1, vs1, r2, vs2))
        
    return layers_p, layers_s

def shoot_analytical(takeoff_angle_deg, source_depth, layers):
    r_src = R_EARTH - source_depth
    # Безопасно определяем скорость источника (если V=0, например S волна в ядре, то 0)
    v_src = next(v1 * (r_src / r1)**(np.log(v2/v1)/np.log(r2/r1) if r1!=r2 and v1>0 and v2>0 else 0) 
                 for (r1, v1, r2, v2) in layers if r1 >= r_src >= r2)
    
    theta_rad = np.radians(takeoff_angle_deg)
    p = (r_src * np.sin(theta_rad)) / v_src if v_src > 0 else 0
    
    path_r, path_delta = [r_src], [0.0]
    current_delta = 0.0
    
    layer_idx = next(i for i, l in enumerate(layers) if l[0] > r_src)
    turn_layer_idx, turn_r = None, None
    reflected_at_boundary = False
    
    # Трассировка ВНИЗ
    for i in range(layer_idx, len(layers)):
        r1, v1, r2, v2 = layers[i]
        if i == layer_idx and r_src < r1: r1, v1 = r_src, v_src
        
        # Если волна вошла в слой где ее скорость <= 0 (S-волна во внешнем ядре) -> отражаем
        if v1 <= 0 or v2 <= 0:
            turn_r, turn_layer_idx, reflected_at_boundary = r1, i - 1, True
            break
            
        zeta = np.log(v2/v1) / np.log(r2/r1)
        arg1, arg2 = p * v1 / r1, p * v2 / r2
        
        if arg1 >= 1.0:
            turn_r, turn_layer_idx, reflected_at_boundary = r1, i - 1, True
            break
        if arg2 >= 1.0:
            turn_r = ( (r1**zeta) / (p * v1) ) ** (1.0 / (zeta - 1.0))
            current_delta += (1.0 / (1.0 - zeta)) * (np.pi/2.0 - np.arcsin(arg1))
            path_r.append(turn_r); path_delta.append(current_delta)
            turn_layer_idx = i
            break
        else:
            current_delta += (1.0 / (1.0 - zeta)) * (np.arcsin(arg2) - np.arcsin(arg1))
            path_r.append(r2); path_delta.append(current_delta)
            
    if turn_layer_idx is None:
        # Для идеально вертикального луча p=0 (takeoff_angle=0)
        if abs(takeoff_angle_deg) < 1e-5:
            turn_layer_idx = len(layers) - 1
            turn_r = layers[-1][2]
            reflected_at_boundary = False
            # Проход напрямую через центр Земли — прыжок азимута ровно на 180°
            current_delta = np.pi 
        else:
            raise ValueError("Луч прошел ядро насквозь")
        
    # Трассировка ВВЕРХ
    for j in range(turn_layer_idx, -1, -1):
        r1_orig, v1_orig, r2_orig, v2_orig = layers[j]
        # Пропускаем слои, где среда не поддерживает тип волны (v=0)
        if v1_orig <= 0 or v2_orig <= 0: continue
            
        zeta = np.log(v2_orig/v1_orig) / np.log(r2_orig/r1_orig)
        r_top, v_top, r_bot, v_bot = r1_orig, v1_orig, r2_orig, v2_orig
        
        # Безопасный расчет аргументов для arcsin даже при p=0
        if j == turn_layer_idx:
            r_bot = turn_r
            if abs(p) < 1e-12:
                arg_bot = 0.0
            else:
                arg_bot = min(p * v_bot / r_bot, 1.0) if reflected_at_boundary else 1.0
        else:
            arg_bot = min(p * v_bot / r_bot, 1.0) if abs(p) > 1e-12 else 0.0
            
        arg_top = min(p * v_top / r_top, 1.0) if abs(p) > 1e-12 else 0.0
        
        current_delta += (1.0 / (1.0 - zeta)) * (np.arcsin(arg_bot) - np.arcsin(arg_top))
        path_r.append(r_top); path_delta.append(current_delta)
        
    return np.degrees(current_delta), path_r, path_delta

def find_all_takeoff_angles(shortest_target, source_depth, layers):
    """ Находит все углы для математически кратчайшей дальности """
    # Запускаем от 0.0 (строго вертикально!), чтобы найти вертикальный луч!
    thetas = np.linspace(0.0, 89.9, 500)
    deltas = []
    for t in thetas:
        try:
            d = shoot_analytical(t, source_depth, layers)[0]
            deltas.append(d)
        except ValueError:
            deltas.append(-999.0)
            
    deltas = np.array(deltas)
    results = []
    
    # Ищем лучи, которые летят коротким путем и длинным путем (360 - dist)
    targets = [shortest_target]
    if shortest_target != 180 and shortest_target != 0:
        targets.append(360.0 - shortest_target)
        
    for target in targets:
        for i in range(len(thetas)-1):
            d1, d2 = deltas[i], deltas[i+1]
            if d1 < -900 or d2 < -900: continue
            
            if (d1 - target) * (d2 - target) <= 0 and abs(d1 - d2) < 20.0:
                func = lambda t: shoot_analytical(t, source_depth, layers)[0] - target
                try:
                    exact_theta = brentq(func, thetas[i], thetas[i+1])
                    results.append((exact_theta, target)) 
                except ValueError:
                    pass
    return results

# ==========================================
# 2. ВИЗУАЛИЗАЦИЯ MANIM
# ==========================================

class SeismicRayTracing(Scene):
    
    def construct(self):
        self.camera.background_color = WHITE
        SCALE = 3.0 / R_EARTH
        
        def pol2cart(r, delta_deg):
            # Жесткое ограничение, чтобы линии из-за погрешности float не вылезали за окружность Земли
            r = min(r, R_EARTH)
            delta_rad = np.radians(delta_deg)
            x = r * SCALE * np.sin(delta_rad)
            y = r * SCALE * np.cos(delta_rad)
            return np.array([x, y, 0])

        # --- ОТРИСОВКА ПЛАНЕТЫ ---
        earth = Circle(radius=R_EARTH * SCALE, color=ORANGE, stroke_width=2)
        earth_fill = Circle(radius=R_EARTH * SCALE, color=ORANGE, fill_opacity=0.1, stroke_width=0)
        
        outer_core = Circle(radius=(R_EARTH - 2889) * SCALE, color=ORANGE, fill_opacity=0.2, stroke_width=2)
        DashedVMobject(outer_core, num_dashes=60)
        inner_core = Circle(radius=(R_EARTH - 5153) * SCALE, color=RED, fill_opacity=0.4, stroke_width=1)
        
        ev_depth = 150.0
        
        source_star = Star(color=RED, fill_opacity=1).scale(0.2).move_to(pol2cart(R_EARTH - ev_depth, 0))

        layers_p, layers_s = build_layers("/home/ksukskss/Projects/NERV/data/IASP91.csv")

        # --- СТАНЦИИ С ОБРАТНЫМ АЗИМУТОМ ---
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
            
            # 1. Сначала ставим маркер станции на поверхность Земли
            st_pos = pol2cart(R_EARTH, target_angle)
            st_marker = Triangle(color=GRAY, fill_opacity=1).scale(0.1).move_to(st_pos)
            
            # Поворачиваем треугольник вершиной наружу
            # (Если Triangle в Manim изначально смотрит вверх, то этот поворот сориентирует его радиально)
            st_marker.rotate(-np.radians(target_angle))
            
            # 2. Считаем вектор направления от центра Земли (0,0) к станции наружу
            direction_vector = normalize(st_marker.get_center())

            # 3. Создаем текст названия
            st_label = Text(st_name, font_size=16, color=GRAY)
            
            # Поворот текста, чтобы он читался комфортно (не вверх ногами)
            st_label.rotate(-np.radians(target_angle))
            st_label.rotate(PI / 2)
            if 180 < target_angle < 360:
                st_label.rotate(PI)

            # 4. ПРИВЯЗКА К МАРКЕРУ:
            # Ставим текст СЛЕДУЮЩИМ ЗА маркером в направлении нашего вектора.
            # buff=0.15 — это размер отступа от вершины треугольника до текста в единицах Manim.
            st_label.next_to(st_marker, direction_vector, buff=0.15)

            # 5. Сборка в группу для корректного отображения слоев (z_index)
            station_group = VGroup(st_marker, st_label)
            station_mobjects.append(station_group)

            # === Расчет лучей P-ВЫСТРЕЛОВ ===
            ray_results_p = find_all_takeoff_angles(shortest_dist, ev_depth, layers_p)
            for angle, traveled_dist in ray_results_p:
                _, path_r, path_delta = shoot_analytical(angle, ev_depth, layers_p)
                
                # ФИЛЬТР ОТ "ПУЧКОВ": отбрасываем лучи, которые не дошли до поверхности
                if abs(path_r[-1] - R_EARTH) > 5.0:
                    continue

                is_short_ray = traveled_dist <= 180
                multiplier = 1 if (is_right_side == is_short_ray) else -1
                points = [pol2cart(r, np.degrees(d) * multiplier) for r, d in zip(path_r, path_delta)]
                
                # Цвет P-волн сделал синим для контраста на белом фоне
                ray_line = VMobject().set_points_as_corners(points).set_color(YELLOW).set_stroke(width=2.5)
                p_rays.append(ray_line)

            # === Расчет лучей S-ВЫСТРЕЛОВ ===
            ray_results_s = find_all_takeoff_angles(shortest_dist, ev_depth, layers_s)
            for angle, traveled_dist in ray_results_s:
                _, path_r, path_delta = shoot_analytical(angle, ev_depth, layers_s)
                
                # ФИЛЬТР ОТ "ПУЧКОВ" для S-волн
                if abs(path_r[-1] - R_EARTH) > 5.0:
                    continue

                is_short_ray = traveled_dist <= 180
                multiplier = 1 if (is_right_side == is_short_ray) else -1
                points = [pol2cart(r, np.degrees(d) * multiplier) for r, d in zip(path_r, path_delta)]

                # 1. Создаем обычную сплошную траекторию луча
                solid_ray_s = VMobject().set_points_as_corners(points)
                
                # 2. Передаем её в DashedVMobject, указав число штрихов через num_dashes
                ray_line_s = DashedVMobject(solid_ray_s, num_dashes=25)
                
                # 3. Задаем цвет и толщину пунктирного луча
                ray_line_s.set_color(BLACK).set_stroke(width=3.0)
                
                s_rays.append(ray_line_s)

        # 1. МГНОВЕННОЕ (СТАТИЧНОЕ) ПОЯВЛЕНИЕ ПЛАНЕТЫ И СТАНЦИЙ
        # self.add() просто кладет объекты на холст с 0-й секунды без анимации
        source_star.set_z_index(10)
        for st in station_mobjects:
            st.set_z_index(10)
        self.add(earth, earth_fill, outer_core, inner_core, source_star, *station_mobjects)
        self.wait(0.5)

        # 2. РАВНОМЕРНАЯ СКОРОСТЬ С РАЗНЫМ ТАЙМИНГОМ
        # rate_func=linear заставляет волну "ползти" равномерно, не ускоряясь в ядре
        # .set_run_time() задает строгую длительность: 2.5 сек для P-волн, 5.0 сек для S-волн
        p_animations = [Create(ray, rate_func=linear).set_run_time(4.0) for ray in p_rays]
        s_animations = [Create(ray, rate_func=linear).set_run_time(4.0) for ray in s_rays]
        
        # Запускаем всё одновременно
        self.play(
            *p_animations,
            *s_animations
        )
            
        self.wait(3)