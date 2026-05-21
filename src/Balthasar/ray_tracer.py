import pandas as pd
import numpy as np

class BalthasarRayTracer:
    def __init__(self, iasp91_csv_path):
        """
        Инициализация лучевого трассировщика под ваш файл IASP91.csv.
        Структура файла: без заголовков, столбцы [depth, radius, vp, vs], шаг ровно 1 км.
        """
        # Читаем ваш файл без заголовков и даем столбцам правильные имена
        df = pd.read_csv(iasp91_csv_path, header=None, names=['depth', 'radius', 'vp', 'vs'])
        
        self.R_E = 6371.0  # Радиус Земли в км
        
        # Вырезаем мантию Земли (до границы с внешним ядром на глубине 2889 км)
        # Это гарантирует, что прямая P-волна не зайдет в ядро
        mantle_df = df[df['depth'] <= 2889.0]
        
        # Загружаем реальные данные шага сетки 1 км напрямую без интерполяции
        self.target_depths = mantle_df['depth'].values
        self.radii = mantle_df['radius'].values
        self.vp = mantle_df['vp'].values
        self.vs = mantle_df['vs'].values

    def compute_leg_to_turn(self, p, r_start, start_idx):
        """
        Аналитическое интегрирование траектории луча (одного плеча) до точки разворота.
        Исправлено: корректно обрабатывает полное внутреннее отражение на границах слоев.
        """
        delta = 0.0
        time = 0.0
        current_r = r_start
        
        for i in range(start_idx, len(self.radii) - 1):
            v = self.vp[i]
            r_bot = self.radii[i+1]
            r_turn = p * v
            
            # --- ИСПРАВЛЕНИЕ ЛОВУШКИ ---
            # Если из-за градиента скорости r_turn превысил или сравнялся с текущим радиусом,
            # это означает ПОЛНОЕ ВНУТРЕННЕЕ ОТРАЖЕНИЕ на границе дискретного слоя.
            if r_turn >= current_r:
                if i == start_idx:
                    return None, None  # Физически невозможный угол выхода из самого источника
                else:
                    return delta, time  # Успешный разворот на границе раздела слоев
            
            # Проверяем, разворачивается ли луч внутри самого 1-км слоя
            if r_turn >= r_bot:
                i1 = np.arcsin(r_turn / current_r)
                i2 = np.pi / 2.0  # Угол падения становится 90 градусов (горизонт)
                
                delta += (i2 - i1)
                time += np.sqrt(current_r**2 - r_turn**2) / v
                return delta, time
            else:
                # Луч беспрепятственно проходит слой насквозь от top (current_r) до bottom (r_bot)
                i1 = np.arcsin(r_turn / current_r)
                i2 = np.arcsin(r_turn / r_bot)
                
                delta += (i2 - i1)
                time += (np.sqrt(current_r**2 - r_turn**2) - np.sqrt(r_bot**2 - r_turn**2)) / v
                current_r = r_bot
                
        return None, None  # Луч пробил мантию насквозь и упал во внешнее ядро

    def shoot_ray(self, p, source_depth):
        r_source = self.R_E - source_depth
        idx_s = int(np.floor(source_depth))
        
        if idx_s >= len(self.radii):
            return None, None, None, None
        
        # Плечо 1: От источника (глубина очага) вниз до точки разворота в мантии
        res_down = self.compute_leg_to_turn(p, r_source, idx_s)
        if res_down[0] is None:
            return None, None, None, None
            
        # Плечо 2: От поверхности Земли до точки разворота (восходящее плечо)
        res_up = self.compute_leg_to_turn(p, self.R_E, 0)
        if res_up[0] is None:
            return None, None, None, None
            
        # Полный путь луча — это сумма нисходящего и восходящего движений
        total_delta = res_down[0] + res_up[0]
        total_time = res_down[1] + res_up[1]
        
        # Вычисление углов выхода и входа по закону Снеллиуса
        v_source = self.vp[idx_s]
        v_surface = self.vp[0]
        
        i_h = np.degrees(np.arcsin(np.clip(p * v_source / r_source, -1.0, 1.0)))
        i_rec = np.degrees(np.arcsin(np.clip(p * v_surface / self.R_E, -1.0, 1.0)))
        
        return total_delta, total_time, i_h, i_rec

    def find_ray_for_station(self, target_delta_rad, source_depth, tolerance=1e-4, max_iter=100):
        idx_s = int(np.floor(source_depth))
        if idx_s >= len(self.radii):
            raise ValueError(f"Глубина источника {source_depth} км выходит за пределы мантии.")
            
        v_source = self.vp[idx_s]
        r_source = self.R_E - source_depth
        
        # Стартовые рамки бисекции
        p_min = 0.0  # Предельно вертикальный луч
        p_max = r_source / v_source - 1e-4  # Предельно горизонтальный луч
        
        for iteration in range(max_iter):
            p_mid = (p_min + p_max) / 2.0
            res = self.shoot_ray(p_mid, source_depth)
            
            if res[0] is None:
                # Если трассировка вернула None, значит мы пробили ядро.
                # Чтобы сделать луч более пологим и удержать в мантии -> увеличиваем p
                p_min = p_mid
                continue
                
            current_delta, current_time, i_h, i_rec = res
            error = current_delta - target_delta_rad
            
            if abs(error) < tolerance:
                return {
                    "ray_parameter": p_mid,
                    "travel_time": current_time,
                    "takeoff_angle": i_h,
                    "incidence_angle": i_rec,
                    "error_deg": np.degrees(error)
                }
            
            # Сейсмическая физика: чем больше параметр p, тем более пологий луч и тем МЕНЬШЕ дистанция
            if error > 0:  
                # Получили дистанцию больше, чем надо -> нужно уменьшить её -> увеличиваем p
                p_min = p_mid
            else:          
                # Получили дистанцию меньше, чем надо -> нужно увеличить её -> уменьшаем p
                p_max = p_mid
                
        raise ValueError(f"[Magi Balthazar Error]: Не сошелся на расстоянии {np.degrees(target_delta_rad):.2f}° за {max_iter} итераций.")