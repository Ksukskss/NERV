#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

const double R_EARTH = 6371.0;
const double PI = std::acos(-1.0);

struct Layer {
    double r1;
    double v1;
    double r2;
    double v2;
    
    // Явный конструктор для безопасной инициализации из Pybind11
    Layer(double r1_val, double v1_val, double r2_val, double v2_val) 
        : r1(r1_val), v1(v1_val), r2(r2_val), v2(v2_val) {}
};

struct RayPath {
    double delta_deg;
    std::vector<double> r;
    std::vector<double> delta;
};

struct RayResult {
    double takeoff_angle;
    double traveled_dist;
    RayPath path;
};

double brentq(std::function<double(double)> f, double a, double b, double tol = 1e-5, int max_iter = 100) {
    double fa = f(a), fb = f(b);
    if (fa * fb >= 0) throw std::runtime_error("Root is not bracketed in brentq!");
    if (std::abs(fa) < std::abs(fb)) { std::swap(a, b); std::swap(fa, fb); }
    
    double c = a, fc = fa, d = 0.0, e = 0.0;
    bool mflag = true;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        if (fb == 0.0 || std::abs(b - a) < tol) return b;
        
        double s;
        if (fa != fc && fb != fc) {
            s = a * fb * fc / ((fa - fb) * (fa - fc)) +
                b * fa * fc / ((fb - fa) * (fb - fc)) +
                c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            s = b - fb * (b - a) / (fb - fa);
        }
        
        double cond1 = (3.0 * a + b) / 4.0;
        bool is_s_between = (s > std::min(cond1, b)) && (s < std::max(cond1, b));
        
        if (!is_s_between ||
            (mflag && std::abs(s - b) >= std::abs(b - c) / 2.0) ||
            (!mflag && std::abs(s - b) >= std::abs(c - d) / 2.0) ||
            (mflag && std::abs(b - c) < tol) ||
            (!mflag && std::abs(c - d) < tol)) {
            s = (a + b) / 2.0; mflag = true;
        } else { mflag = false; }
        
        double fs = f(s);
        d = c; c = b; fc = fb;
        
        if (fa * fs < 0) { b = s; fb = fs; } 
        else { a = s; fa = fs; }
        
        if (std::abs(fa) < std::abs(fb)) { std::swap(a, b); std::swap(fa, fb); }
    }
    return b;
}

RayPath shoot_analytical(double takeoff_angle_deg, double source_depth, const std::vector<Layer>& layers, bool return_path = true) {
    if (layers.empty()) throw std::runtime_error("Earth model layers are empty!");
    
    double r_src = R_EARTH - source_depth;
    
    // Защита от выхода за границы модели Земли
    if (r_src > layers.front().r1) r_src = layers.front().r1;
    if (r_src < layers.back().r2) r_src = layers.back().r2;

    double v_src = 0.0;
    int src_layer_idx = -1;
    
    for (size_t i = 0; i < layers.size(); ++i) {
        // Добавлен допуск 1e-7 для избежания проблем с точностью float
        if (layers[i].r1 + 1e-7 >= r_src && r_src >= layers[i].r2 - 1e-7) {
            src_layer_idx = static_cast<int>(i);
            if (layers[i].v1 > 0 && layers[i].v2 > 0 && std::abs(layers[i].r1 - layers[i].r2) > 1e-9) {
                double zeta = std::log(layers[i].v2 / layers[i].v1) / std::log(layers[i].r2 / layers[i].r1);
                v_src = layers[i].v1 * std::pow(r_src / layers[i].r1, zeta);
            }
            break;
        }
    }
    
    if (src_layer_idx == -1) {
        throw std::runtime_error("Source depth does not match any layer.");
    }
    
    double theta_rad = takeoff_angle_deg * PI / 180.0;
    double p = (v_src > 0.0) ? (r_src * std::sin(theta_rad)) / v_src : 0.0;
    
    RayPath result;
    double current_delta = 0.0;
    
    if (return_path) {
        result.r.push_back(r_src);
        result.delta.push_back(0.0);
    }
    
    int turn_layer_idx = -1;
    double turn_r = -1.0;
    bool reflected_at_boundary = false;
    
    for (int i = src_layer_idx; i < static_cast<int>(layers.size()); ++i) {
        double r1 = layers[i].r1;
        double v1 = layers[i].v1;
        double r2 = layers[i].r2;
        double v2 = layers[i].v2;
        
        if (i == src_layer_idx && r_src < r1) { r1 = r_src; v1 = v_src; }

        // Skip a source layer that has collapsed to zero thickness
        // (happens when the source sits exactly on a layer boundary).
        if (std::abs(r1 - r2) < 1e-7) { continue; }

        if (v1 <= 0.0 || v2 <= 0.0) {
            turn_r = r1; turn_layer_idx = i - 1; reflected_at_boundary = true; break;
        }
        
        double zeta = std::log(v2 / v1) / std::log(r2 / r1);
        double denom = 1.0 - zeta;
        // Защита от деления на ноль
        if (std::abs(denom) < 1e-9) denom = (denom >= 0) ? 1e-9 : -1e-9;
        
        double arg1 = p * v1 / r1;
        double arg2 = p * v2 / r2;
        
        if (arg1 >= 1.0) {
            turn_r = r1; turn_layer_idx = i - 1; reflected_at_boundary = true; break;
        }
        
        if (arg2 >= 1.0) {
            turn_r = std::pow( std::pow(r1, zeta) / (p * v1), 1.0 / -denom );
            current_delta += (1.0 / denom) * (PI / 2.0 - std::asin(std::min(1.0, arg1)));
            if (return_path) { result.r.push_back(turn_r); result.delta.push_back(current_delta); }
            turn_layer_idx = i; break;
        } else {
            current_delta += (1.0 / denom) * (std::asin(std::min(1.0, arg2)) - std::asin(std::min(1.0, arg1)));
            if (return_path) { result.r.push_back(r2); result.delta.push_back(current_delta); }
        }
    }
    
    if (turn_layer_idx == -1) {
        if (std::abs(takeoff_angle_deg) < 1e-5) {
            turn_layer_idx = static_cast<int>(layers.size()) - 1; 
            turn_r = layers.back().r2; 
            reflected_at_boundary = false; 
            current_delta = PI; 
        } else { 
            throw std::runtime_error("Ray penetrated the core completely."); 
        }
    }
    
    for (int j = turn_layer_idx; j >= 0; --j) {
        double r1_orig = layers[j].r1, v1_orig = layers[j].v1, r2_orig = layers[j].r2, v2_orig = layers[j].v2;
        if (v1_orig <= 0.0 || v2_orig <= 0.0) continue;
        
        double zeta = std::log(v2_orig / v1_orig) / std::log(r2_orig / r1_orig);
        double denom = 1.0 - zeta;
        if (std::abs(denom) < 1e-9) denom = (denom >= 0) ? 1e-9 : -1e-9;
        
        double arg_bot = 0.0, arg_top = 0.0;
        
        if (j == turn_layer_idx) {
            if (std::abs(p) > 1e-12) arg_bot = reflected_at_boundary ? std::min(p * v2_orig / turn_r, 1.0) : 1.0;
        } else {
            if (std::abs(p) > 1e-12) arg_bot = std::min(p * v2_orig / r2_orig, 1.0);
        }
        
        if (std::abs(p) > 1e-12) arg_top = std::min(p * v1_orig / r1_orig, 1.0);
        
        current_delta += (1.0 / denom) * (std::asin(arg_bot) - std::asin(arg_top));
        if (return_path) { result.r.push_back(r1_orig); result.delta.push_back(current_delta); }
    }
    
    result.delta_deg = current_delta * 180.0 / PI;
    return result;
}

std::vector<RayResult> find_all_takeoff_angles(double shortest_target, double source_depth, const std::vector<Layer>& layers) {
    int num_steps = 500;
    std::vector<double> thetas(num_steps);
    std::vector<double> deltas(num_steps);
    
    for (int i = 0; i < num_steps; ++i) {
        thetas[i] = 0.0 + i * (89.9 - 0.0) / (num_steps - 1);
        try { 
            deltas[i] = shoot_analytical(thetas[i], source_depth, layers, false).delta_deg; 
        } 
        catch (...) { 
            deltas[i] = -999.0; 
        }
    }
    
    std::vector<double> targets = {shortest_target};
    if (shortest_target != 180.0 && shortest_target != 0.0) targets.push_back(360.0 - shortest_target);
    
    std::vector<RayResult> results;
    for (double target : targets) {
        for (int i = 0; i < num_steps - 1; ++i) {
            double d1 = deltas[i], d2 = deltas[i+1];
            if (d1 < -900 || d2 < -900) continue;
            
            if ((d1 - target) * (d2 - target) <= 0.0 && std::abs(d1 - d2) < 20.0) {
                auto func = [&](double t) { return shoot_analytical(t, source_depth, layers, false).delta_deg - target; };
                try {
                    double exact_theta = brentq(func, thetas[i], thetas[i+1]);
                    RayPath full_path = shoot_analytical(exact_theta, source_depth, layers, true);
                    results.push_back({exact_theta, target, full_path});
                } catch (const std::exception& e) {
                    std::cerr << "C++ Error at target " << target << " deg: " << e.what() << "\n";
                }
            }
        }
    }
    
    // Печатаем прямо в консоль для отладки
    if (!results.empty()) {
        std::cout << "C++: Target " << shortest_target << " deg: Found " << results.size() << " valid rays.\n";
    }
    
    return results;
}

PYBIND11_MODULE(ray_tracer_cpp, m) {
    m.doc() = "C++ plugin for fast seismic ray tracing";

    py::class_<Layer>(m, "Layer")
        .def(py::init<double, double, double, double>())
        .def_readwrite("r1", &Layer::r1)
        .def_readwrite("v1", &Layer::v1)
        .def_readwrite("r2", &Layer::r2)
        .def_readwrite("v2", &Layer::v2);

    py::class_<RayPath>(m, "RayPath")
        .def_readonly("delta_deg", &RayPath::delta_deg)
        .def_readonly("r", &RayPath::r)
        .def_readonly("delta", &RayPath::delta);

    py::class_<RayResult>(m, "RayResult")
        .def_readonly("takeoff_angle", &RayResult::takeoff_angle)
        .def_readonly("traveled_dist", &RayResult::traveled_dist)
        .def_readonly("path", &RayResult::path);

    m.def("find_all_takeoff_angles", &find_all_takeoff_angles, "Find takeoff angles and return exact paths");
}