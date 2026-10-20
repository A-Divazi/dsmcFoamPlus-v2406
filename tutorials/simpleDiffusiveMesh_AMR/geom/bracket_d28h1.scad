// Change this to false to remove the helper geometry
debug = false;

// Global resolution
$fs = 0.1;  // Don't generate smaller facets than 0.1 mm
$fa = 0.1;    // Don't generate larger angles than 5 degrees
$fn=100;


module circle_ring(radius, thickness) {
    rotate_extrude(angle=360, convexity = 0){
                        translate([radius, 0, 0])
                        square(size = [thickness, thickness], center = true);
        }
}

scale(1e-3) translate([0.5, 0, 0]) rotate([0, 90, 0]) circle_ring(13.5, 1); //see fig.1 in Gumbel, AGP 2001: https://doi.org/10.5194/angeo-19-563-2001