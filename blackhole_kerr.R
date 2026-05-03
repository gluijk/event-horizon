# Black hole basic ray tracing simulation
# Improved version with: normalised FOV, antialiasing and approximated Kerr generalisation
# www.overfitting.net
# https://www.overfitting.net/2026/04/simulando-agujeros-negros.html

# Other resources:
# History of first black hole simulation:
#   https://astrobitos.org/2019/05/02/la-historia-detras-la-primera-imagen-simulada-de-un-agujero-negro/
#   https://www.cnrs.fr/en/press/first-ever-image-black-hole-cnrs-researcher-had-simulated-it-early-1979
# Realtime Javascript black hole raytracer:
#   https://adriwin06.github.io/black-hole/


library(tiff)
library(png)
library(Rcpp)


# How to Tweak the Results:

# Resolution: Adjust the width and height parameters. Higher resolutions take linearly
# more time to render

# Camera Angle: Modify cam_elev. 0.0 will give you a dead-on edge view,
# while 1.5 will give you a top-down view where the gravitational lensing is
# less pronounced but the orbit is clearer


# 1. Define the C++ raytracer (courtesy of Gemini Pro)

# BASIC RAYTRACER + KERR BLACK HOLE APPROX + AA: render_bh_cpp_kerr()
# ChatGPT verdict: this is a visually plausible renderer, not a physically accurate Kerr ray tracer
cpp_code <- "
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <random>

using namespace Rcpp;

// Helper: Vector cross product
inline void cross(const double a[3], const double b[3], double res[3]) {
    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
}

// Helper: Vector dot product
inline double dot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Helper: Normalize vector
inline void normalize(double v[3]) {
    double len = std::sqrt(dot(v, v));
    if (len > 0) { v[0]/=len; v[1]/=len; v[2]/=len; }
}

// Helper: Signum function
inline double sgn(double val) {
    return (double(0.0) < val) - (val < double(0.0));
}

// [[Rcpp::export]]
NumericVector render_bh_cpp_kerr(
        int width=1920, int height=1080,                                    // output resolution in pixels
        double cam_dist=50, double cam_elev=0.075, double FOV_scale = 0.6,  // camera parameters
        double a_star=0.0, double M=1.0, double r_out_accretion = 20.0,     // black hole parameters
        int glow=1, int rings=1, int AA=1) {                                // plot parameters

    // cam_dist: camera distance in M units
    // cam_elev: camera elevation in radians
    // FOV_scale: viewport's physical full height relative to a focal length of 1
    //            VFOV = 2 * atan(0.6 / 2) ~33.4º / HFOV = 2 * atan(0.6 * aspect_ratio / 2)
    // a_star: user input spin parameter is assumed to be normalized a* in the {-1,+1} dimensionless range
    //         while the formulas are referred to a in the {−M,+M} units range
    //         a* = a / M -> a = a* * M
    // M: black hole mass in length units
    // r_out_accretion: size of outer radius for accretion disk in length units
    // glow: boolean to choose colour palette
    // rings: boolean to add rings on accretion disk
    // AA: number of antialiasing pixels (AA=1 means no antialiasing)

    double a = a_star * M;  // scaling to give a units of M (length) which all used formulas expect

    NumericVector img(width * height * 3);
    
    // Camera setup
    double cx = 0.0;
    double cy = -cam_dist * std::cos(cam_elev);
    double cz = cam_dist * std::sin(cam_elev);
    
    double dt = 0.1;           

    // Antialiasing random generator
    std::mt19937 rng(42); // seed
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Calculate Event horizon
    double r_H = M + std::sqrt(std::max(0.0, M*M - a*a));  // independent of the sign of a

    // Calculate ISCO for Kerr (formula by Bardeen et al, 1972)
    // if a=0 (Schwarzschild) -> Z1=3, Z2=3, r_ISCO=6
    double Z1 = 1.0 + std::pow(1.0 - a*a/(M*M), 1.0/3.0) * (std::pow(1.0 + a/M, 1.0/3.0) + std::pow(1.0 - a/M, 1.0/3.0));
    double Z2 = std::sqrt(3.0 * a*a/(M*M) + Z1*Z1);
    double r_ISCO = M * (3.0 + Z2 - sgn(a) * std::sqrt((3.0 - Z1)*(3.0 + Z1 + 2.0*Z2)));  // sgn(a) accounts for prograde/retrograde

    // Display main black hole parameters and radius
    const char* label = (a == 0) 
        ? \"Schwarzschild black hole: \" 
        : \"Kerr black hole: \";
    Rcpp::Rcout << label
                << \"M=\" << M
                << \", a*=\" << a_star
                << \" -> a=\" << a
                << \", r_H=\" << r_H
                << \", r_ISCO=\" << r_ISCO
                << std::endl;

    double aspect_ratio = (double)width / (double)height;
    double cdir[3] = {0.0, -cy, -cz}; 
    normalize(cdir);
    double right[3] = {1.0, 0.0, 0.0};
    double up[3];
    cross(right, cdir, up);

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {

            // Pixel colour accumulator
            double col_acc[3] = {0.0, 0.0, 0.0};
            
            // ANTIALIASING LOOP: send N rays per pixel
            for (int s = 0; s < AA; ++s) {
                
                // Add random offset between 0 and 1 (or 0.5 if AA=1)
                double offset_u = (AA == 1) ? 0.5 : dist(rng);
                double offset_v = (AA == 1) ? 0.5 : dist(rng);
                
                // Map pixel to screen space coordinates
                double u_scr = (double(i + offset_u) / width - 0.5) * FOV_scale * aspect_ratio;
                // Invert Y so the top of the image renders correctly
                double v_scr = -(double(j + offset_v) / height - 0.5) * FOV_scale;

                double pos[3] = {cx, cy, cz};
                double vel[3] = {
                    cdir[0] + u_scr*right[0] + v_scr*up[0],
                    cdir[1] + u_scr*right[1] + v_scr*up[1],
                    cdir[2] + u_scr*right[2] + v_scr*up[2]
                };
                normalize(vel);
                
                double color[3] = {0.0, 0.0, 0.0};
                double z_prev = pos[2];
                
                for (int step = 0; step < 1000; ++step) {
                    double r2 = dot(pos, pos);
                    double r = std::sqrt(r2);
                    
                    if (r < r_H * 1.01) break;  // early stop 1% safety margin near event horizon (visually indistinguishable)
                    if (r > cam_dist * 2.0) break;  // 2.0 factor to give up on a photon that is flying away into deep space
                    
                    if (pos[2] * z_prev <= 0 && step > 0) {  // sign change in z: the ray has passed through z = 0
                        double t = std::abs(z_prev) / (std::abs(z_prev) + std::abs(pos[2]) + 1e-8);
                        // NOTE: hit_x, hit_y and hit_r are coming from the crossing point between
                        // the photon trajectory (ray segment between two integration steps) and the equatorial plane z=0
                        double hit_x = pos[0] - vel[0]*dt*(1.0-t);
                        double hit_y = pos[1] - vel[1]*dt*(1.0-t);
                        double hit_r = std::sqrt(hit_x*hit_x + hit_y*hit_y);
                        
                        if (hit_r > r_ISCO && hit_r < r_out_accretion) {
                            // Keplerian velocity: Omega = 1 / (r^1.5 + a) naturally handles the sign of a
                            double omega = 1.0 / (std::pow(hit_r, 1.5)/std::sqrt(M) + a);
                            double vk = omega * hit_r;
                            
                            double V_disk[3] = {-vk * hit_y / hit_r, vk * hit_x / hit_r, 0.0};
                            double n[3] = {-vel[0], -vel[1], -vel[2]};
                            
                            double V_dot_n = dot(V_disk, n);
                            double g = std::sqrt(1.0 - vk*vk) / (1.0 - V_dot_n);
                            
                            // Redshift approx
                            double delta = hit_r*hit_r - 2.0*M*hit_r + a*a;
                            double sigma = hit_r*hit_r; 
                            double red_grav = std::sqrt(delta / sigma);
                            g *= red_grav;
                            
                            if (glow) {  // orangish colours
                                double base_I = 150.0 / (hit_r * hit_r);  // 1/r2 law
                                double I_obs = base_I * std::pow(g, 4.0);
                                color[0] = std::min(1.0, I_obs * 1.0 * std::pow(g, 1.5)); // Red
                                color[1] = std::min(1.0, I_obs * 0.6 * std::pow(g, 3.0)); // Green
                                color[2] = std::min(1.0, I_obs * 0.2 * std::pow(g, 4.5)); // Blue
                            } else {  // bluish/orangish colours
                                double base_I = 200.0 / (std::pow(hit_r, 3.0));  // 1/r3 law
                                double I_obs = base_I * std::pow(g, 4.0);
                                color[0] = std::pow(std::min(1.0, I_obs * std::pow(g, 1.2)), 0.5);
                                color[1] = std::pow(std::min(1.0, I_obs * std::pow(g, 2.5)), 0.5);
                                color[2] = std::pow(std::min(1.0, I_obs * std::pow(g, 4.0)), 0.5);
                            }
    
                            // Optional: add grid rings for depth texture
                            if (rings && std::fmod(hit_r, 1) < 0.1) {  // frequency, pulse width
                                color[0] *= 0.3; color[1] *= 0.3; color[2] *= 0.3;
                            }
                            break;  // stop tracking photon after it hits the opaque disk
                        }
                    }
                    
                    z_prev = pos[2];
                    
                    // --- GEODESIC INTEGRATION ---
                    // Schwarzschild-ish approximation -> photon ring doesn't adapt to spin (a)
                    // so the unstable photon orbit is still effectively at r_ph~3
                    double L[3];
                    cross(pos, vel, L);
                    double L2 = dot(L, L);
                    double a_schwarz = -3.0 * M * L2 / (r2 * r2 * r);
                    
                    // Drag magnitude depends on a, rotation around Z
                    double drag_mag = 2.0 * a * M / (r2 * r);
                    double drag_vec[3] = {-pos[1], pos[0], 0.0}; 
                    
                    vel[0] += (a_schwarz * pos[0] + drag_mag * drag_vec[0]) * dt;
                    vel[1] += (a_schwarz * pos[1] + drag_mag * drag_vec[1]) * dt;
                    vel[2] += (a_schwarz * pos[2]) * dt;
                    normalize(vel);
    
                    pos[0] += vel[0] * dt;
                    pos[1] += vel[1] * dt;
                    pos[2] += vel[2] * dt;
                }

                // Add ray to pixel accumulators
                col_acc[0] += color[0];
                col_acc[1] += color[1];
                col_acc[2] += color[2];
            }

            // Map 1D C++ array to match R's 3D array layout for easy rendering
            img[j + height * i + height * width * 0] = col_acc[0] / AA;
            img[j + height * i + height * width * 1] = col_acc[1] / AA;
            img[j + height * i + height * width * 2] = col_acc[2] / AA;
        }
    }
    
    // Set R array dimensions [height, width, channels]
    img.attr(\"dim\") = IntegerVector::create(height, width, 3);
    return img;
}
"


# 2. Compile the C++ function
cat("Compiling C++ raytracer...\n")
sourceCpp(code = cpp_code)


# Quick default black hole test
img_data <- render_bh_cpp_kerr()
writeTIFF(img_data, "blackhole_default.tif", bits.per.sample = 16)



##############################
# VARYING SPIN ANIMATION


# 3. Setup Camera
OVERSAMPLING=1
width <- 1920*OVERSAMPLING
height <- 1080*OVERSAMPLING
cam_dist <- 30
cam_elev <- 0.075  # angle above the accretion disk (radians)


# 4. Build animation frames
frame=1
for (a in seq(from=-0.98, to=0.98, by=0.02)) {
    name=sprintf("blackhole_%05d.png", frame)
    cat(paste0(name, ": ", sprintf("rendering %dx%d image with a=%f...\n", width, height, a)))
    img_data <- render_bh_cpp_kerr(width, height, cam_dist, cam_elev, FOV_scale=1.2, a_star=a, glow=0)
    writePNG(img_data, name)
    frame=frame+1
}


# Build forward and backward animated 25fps GIF
# magick -delay 4 blackhole_*.png -resize 512x ( -clone 1--2 -reverse ) -loop 0 blackhole.gif



##############################
# BLACK HOLE RADII

library(Cairo)


# --- Black hole radii functions ---
r_g <- function(M, a) { rep(M, length(a)) }  # gravitational radius (r_g)
r_s <- function(M, a) { rep(2 * M, length(a)) }  # Schwarzschild radius (r_s)
r_H <- function(M, a) { M * (1 + sqrt(1 - a^2)) }  # event horizon radius (r_H)
r_ph <- function(M, a) { 2 * M * (1 + cos((2/3) * acos(-a))) }  # photon orbit radius (r_ph)
r_ISCO <- function(M, a) {  # ISCO radius
    Z1 <- 1 + (1 - a^2)^(1/3) * ((1 + a)^(1/3) + (1 - a)^(1/3))
    Z2 <- sqrt(3 * a^2 + Z1^2)
    M * (3 + Z2 - sign(a) * sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)))
}


# --- Calculate radii ---
# Parameters (M, a)
M <- 1
a_vals <- seq(-1, 1, length.out = 801)

rg_vals   <- r_g(M, a_vals)
rH_vals   <- r_H(M, a_vals)
rph_vals  <- r_ph(M, a_vals)
rISCO_vals<- r_ISCO(M, a_vals)
rs_vals   <- r_s(M, a_vals)


# --- Plot radii ---
CairoPNG("blackhole_radii.png", width=512, height=550)
    xlim=c(min(a_vals), max(a_vals))
    ylim= c(0, max(rg_vals, rs_vals, rH_vals, rph_vals, rISCO_vals))
    plot(NA, xlim=xlim, ylim =ylim,
         xlab = "a (spin)", ylab = "Radius / r_g",
         main = "Black hole radii (normalized vs r_g = GM/c^2)",
         yaxt = "n")
    axis(2, at = 0:max(ylim))  # y-axis ticks at every integer
    abline(h=c(3,6,ylim), v=c(0,xlim), col='lightgray')
    
    lines(a_vals, rg_vals, lwd = 2, lty = 3)
    lines(a_vals, rs_vals, lwd = 2, lty = 2)
    lines(a_vals, rH_vals, col = "red", lwd = 2)
    lines(a_vals, rph_vals, col = "darkgreen", lwd = 2)
    lines(a_vals, rISCO_vals, col = "blue", lwd = 2)
    
    legend("topright",
           legend = c("r_ISCO", "Photon orbit (r_ph)", "Event horizon (r_H)",
                      "Schwarzschild (r_s)", "Gravitational (r_g)"),
           col = c("blue", "darkgreen", "red", "black", "black"),
           lty = c(1,1,1,2,3), lwd = 2)
dev.off()
