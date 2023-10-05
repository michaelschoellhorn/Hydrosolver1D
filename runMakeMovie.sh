ffmpeg -framerate 30 -pattern_type glob -i "Plots/q1gauss/*.png" -c:v libx264 -pix_fmt yuv420p gauss_rho.mp4
ffmpeg -framerate 30 -pattern_type glob -i "Plots/q2gauss/*.png" -c:v libx264 -pix_fmt yuv420p gauss_u.mp4
ffmpeg -framerate 30 -pattern_type glob -i "Plots/q1/*.png" -c:v libx264 -pix_fmt yuv420p jump_no_visc_rho.mp4
ffmpeg -framerate 30 -pattern_type glob -i "Plots/q2/*.png" -c:v libx264 -pix_fmt yuv420p jump_no_visc_u.mp4
ffmpeg -framerate 30 -pattern_type glob -i "Plots/q3/*.png" -c:v libx264 -pix_fmt yuv420p jump_no_visc_eps.mp4
ffmpeg -framerate 30 -pattern_type glob -i "Plots/q1visc/*.png" -c:v libx264 -pix_fmt yuv420p jump_visc_rho.mp4
ffmpeg -framerate 30 -pattern_type glob -i "Plots/q2visc/*.png" -c:v libx264 -pix_fmt yuv420p jump_visc_u.mp4
ffmpeg -framerate 30 -pattern_type glob -i "Plots/q3visc/*.png" -c:v libx264 -pix_fmt yuv420p jump_visc_eps.mp4