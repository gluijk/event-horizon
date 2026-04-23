# Black hole basic ray tracing simulation powered by Gemini Pro
# www.overfitting.net
# https://www.overfitting.net/

# Other resources:
# History of first black hole simulation:
#   https://astrobitos.org/2019/05/02/la-historia-detras-la-primera-imagen-simulada-de-un-agujero-negro/
#   https://www.cnrs.fr/en/press/first-ever-image-black-hole-cnrs-researcher-had-simulated-it-early-1979
# Realtime Javascript black hole raytracer:
#   https://adriwin06.github.io/black-hole/


library(tiff)
library(Rcpp)


# How to Tweak the Results:

# Resolution: Adjust the width and height parameters. Higher resolutions take linearly
# more time to render

# Camera Angle: Modify cam_elevation. 0.0 will give you a dead-on edge view,
# while 1.5 will give you a top-down view where the gravitational lensing is
# less pronounced but the orbit is clearer

# Physics: In the C++ code, change r_isco to 3.0 to simulate the photon orbit
# closer to the event horizon, mimicking the physics of a rapidly spinning (Kerr)
# black hole rather than a static (Schwarzschild) one (r_isco=6)
# ISCO=Innermost Stable Circular Orbit


# 1. Define the C++ raytracer
cpp_code <- "
#include <Rcpp.h>
#include <cmath>
#include <algorithm>

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

// [[Rcpp::export]]
NumericVector render_bh_cpp(int width, int height, double cam_dist, double cam_angle, double r_isco) {
    // Array to store RGB values: dims = [height, width, 3]
    NumericVector img(width * height * 3);
    
    // Camera setup (cam_angle in radians. 0 = edge-on, M_PI/2 = top-down)
    double cy = -cam_dist * std::cos(cam_angle);
    double cz = cam_dist * std::sin(cam_angle);
    double cx = 0.0;
    
    // Ray tracing parameters
    double dt = 0.05;         // Integration step size
    double rs = 2.0;          // Schwarzschild radius (M=1)
    // double r_isco = 6.0;      // Innermost Stable Circular Orbit
    double r_out = 20.0;      // Outer edge of the accretion disk
    double fov_scale = 28.0;  // Field of view scalar
    
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            
            // Map pixel to screen space coordinates
            double aspect_ratio = (double)width / (double)height;
            double u = (double(i) / width - 0.5) * fov_scale * aspect_ratio;
            // Invert Y so the top of the image renders correctly
            double v_scr = -(double(j) / height - 0.5) * fov_scale;
            
            // Camera vectors
            double cdir[3] = {0.0, -cy, -cz}; // Pointing to origin
            normalize(cdir);
            double right[3] = {1.0, 0.0, 0.0};
            double up[3];
            cross(right, cdir, up);
            
            double pos[3] = {cx, cy, cz};
            
            // Initial photon velocity (towards screen pixel)
            double vel[3] = {
                cdir[0] + u*right[0]/cam_dist + v_scr*up[0]/cam_dist,
                cdir[1] + u*right[1]/cam_dist + v_scr*up[1]/cam_dist,
                cdir[2] + u*right[2]/cam_dist + v_scr*up[2]/cam_dist
            };
            normalize(vel);
            
            double color[3] = {0.0, 0.0, 0.0}; // Default background: black
            double z_prev = pos[2];
            
            // Ray integration loop
            for (int step = 0; step < 1500; ++step) {
                double r2 = dot(pos, pos);
                double r = std::sqrt(r2);
                
                // Absorbed by black hole
                if (r < rs) break; 
                
                // Escaped gravitational pull
                if (r > cam_dist * 1.5) break; 
                
                // Check if photon crossed the equatorial plane (z = 0)
                if (pos[2] * z_prev <= 0 && step > 0) {
                    
                    // Linear interpolation to exact z=0 crossing point
                    double t = std::abs(z_prev) / (std::abs(z_prev) + std::abs(pos[2]) + 1e-8);
                    double cross_x = pos[0] - vel[0]*dt*(1.0-t);
                    double cross_y = pos[1] - vel[1]*dt*(1.0-t);
                    double cross_r = std::sqrt(cross_x*cross_x + cross_y*cross_y);
                    
                    if (cross_r > r_isco && cross_r < r_out) {
                        // HIT ACCRETION DISK
                        
                        // 1. Relativistic Kinematics
                        // Keplerian velocity: v = sqrt(M/r). Since Rs=2M, M=1. v = sqrt(1/r)
                        double vk = std::sqrt(1.0 / cross_r);
                        
                        // Disk velocity vector (counter-clockwise orbit)
                        double V_disk[3] = {-vk * cross_y / cross_r, vk * cross_x / cross_r, 0.0};
                        
                        // Photon direction is opposite to its velocity
                        double n[3] = {-vel[0], -vel[1], -vel[2]};
                        normalize(n);
                        
                        // 2. Doppler Factor Calculation
                        // g = sqrt(1 - v^2) / (1 - V_disk dot n)
                        double v2 = vk * vk;
                        double V_dot_n = dot(V_disk, n);
                        double g = std::sqrt(1.0 - v2) / (1.0 - V_dot_n);
                        
                        // 3. Emission & Observed Intensity
                        // Intrinsic brightness roughly follows a temperature gradient T ~ r^-0.75 -> I ~ r^-3
                        double base_I = 150.0 / (cross_r * cross_r); 
                        double I_obs = base_I * std::pow(g, 4); // Doppler beaming effect
                        
                        // 4. Color Shifting (Blueshift vs Redshift)
                        double r_col = 1.0 * std::pow(g, 1.5); // Red shifts slower
                        double g_col = 0.6 * std::pow(g, 3.0);
                        double b_col = 0.2 * std::pow(g, 4.5); // Blue spikes on approach
                        
                        color[0] = std::min(1.0, I_obs * r_col);
                        color[1] = std::min(1.0, I_obs * g_col);
                        color[2] = std::min(1.0, I_obs * b_col);
                        
                        // Optional: Add grid rings for depth texture


                        // METHOD 1: just add rings
                        if (std::fmod(cross_r, 1) < 0.1) {  // frequency, pulse width
                            color[0] *= 0.4; color[1] *= 0.4; color[2] *= 0.4;
                        }
                        // To make rings further apart: Increase the 1.5 -> 1 (frequency)
                        // To make rings thinner: Decrease the 0.15 -> 0.1 (pulse width)


                        // To make rings harder or softer: METHOD 1 is a hard cut
                        // If you replaced the if statement with a smooth sin or cos function,
                        // you would get soft, wave-like ripples instead of sharp lines


                        // METHOD 2: rippled rings
                        // 1. Frequency: Multiply cross_r by a higher number to make rings closer (e.g., 4.0)
                        // 2. Softness/Thinness: Use pow() on a shifted sine wave
                        // double ring_pattern = std::sin(cross_r * 4.0); 
                        
                        // Shift sine from [-1, 1] to [0, 1]
                        // ring_pattern = (ring_pattern + 1.0) / 2.0;
                        
                        // Raising to a high power (e.g., 8.0) makes the peak very narrow (thinner)
                        // but keeps the transition smooth (softer)
                        // double darkness_factor = 0.5 + 0.5 * std::pow(ring_pattern, 8.0);
                        
                        // color[0] *= darkness_factor;
                        // color[1] *= darkness_factor;
                        // color[2] *= darkness_factor;
                        

                        break; // Stop tracking photon after it hits the opaque disk
                    }
                }
                
                z_prev = pos[2];
                
                // Calculate photon spatial acceleration (Pseudo-Newtonian GR approximation)
                double L[3];
                cross(pos, vel, L);
                double L2 = dot(L, L);
                double a_coeff = -3.0 * L2 / (r2 * r2 * r);
                
                // Update velocity and position
                vel[0] += a_coeff * pos[0] * dt;
                vel[1] += a_coeff * pos[1] * dt;
                vel[2] += a_coeff * pos[2] * dt;
                normalize(vel); // Force speed of light c=1
                
                pos[0] += vel[0] * dt;
                pos[1] += vel[1] * dt;
                pos[2] += vel[2] * dt;
            }
            
            // Map 1D C++ array to match R's 3D array layout for easy rendering
            img[j + height * i + height * width * 0] = color[0];
            img[j + height * i + height * width * 1] = color[1];
            img[j + height * i + height * width * 2] = color[2];
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


# 3. Setup Camera and Render
OVERSAMPLING=2
width <- 1920*OVERSAMPLING
height <- 1080*OVERSAMPLING
cam_distance <- 50 # 20.0
cam_elevation <- 0.15/2  # angle above the accretion disk (radians). Try 0.4 for a higher view!
r_isco=5.5  # to simulate a photon orbit closer to the Event Horizon set a lower ISCO

cat(sprintf("Rendering %dx%d image...\n", width, height))
system.time({
    # Call the compiled C++ function
    img_data <- render_bh_cpp(width, height, cam_distance, cam_elevation, r_isco)
})


# 4. Save and Display the Image
writeTIFF(img_data, "gemini_r_isco5.5_HQ.tif", bits.per.sample = 16)

# Create a blank plot area, removing margins for a clean image frame
par(mar = c(0, 0, 0, 0), bg = "black")
plot(0, 0, type = 'n', xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = '', ylab = '')
# Draw the raster matrix
rasterImage(as.raster(img_data), 0, 0, 1, 1, interpolate = TRUE)
cat("Done! Notice the bright, blue-shifted left side caused by the relativistic Doppler effect.\n")

