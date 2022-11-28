use image::{ImageBuffer, Rgb};
use std::cmp;

//types
type Color = (f64, f64, f64);
type Point = (f64, f64, f64);
type Normal = (f64, f64, f64);
type Dir = (f64, f64, f64);

type Vector = (f64, f64, f64);


//classes
#[derive(Clone)]
pub struct Plane {
    pub point: Point,
    pub normal: Normal,
    pub color: Color,
    pub diffuse: f64,
    pub reflect: f64,
}

#[derive(Clone)]
pub struct Sphere {
        pub center: Point,
        pub radius: f64,
        pub color: Color,
        pub k_diffuse: f64,
        pub k_reflect: f64,
}

impl Sphere {
    fn getNormal(&self, surfacePoint: Point) -> Normal{
        let dist_from_center = (surfacePoint.0 - self.center.0, surfacePoint.1 - self.center.1, surfacePoint.2 - self.center.2);
        let norm: Normal = (dist_from_center.0 / self.radius, dist_from_center.1 / self.radius, dist_from_center.2 / self.radius);
        return norm;
    }
}


#[derive(Clone)]
pub struct Ray {
    pub point: Point,
    pub dir: Dir,
}



//funcs
fn convert_unit_normal(normal: Normal) -> Normal{
    let mag = (normal.0*normal.0 + normal.1*normal.1 + normal.2*normal.2).sqrt();
    return (normal.0/mag, normal.1/mag, normal.2/mag);
}

fn convert_unit_dir(dir: Dir) -> Dir{
    let mag = (dir.0*dir.0 + dir.1*dir.1 + dir.2*dir.2).sqrt();
    return (dir.0/mag, dir.1/mag, dir.2/mag);
}





fn find_intersection_plane(ray: Ray, plane: Plane, epsilon: f64) -> f64{
    let normalDotRayDir = plane.normal.0*ray.dir.0 + plane.normal.1*ray.dir.1 + plane.normal.2*ray.dir.2;

    if normalDotRayDir.abs() < epsilon {
	return 0.0;
    }

    let normal_unit_vec = convert_unit_normal(plane.normal);
    let distances_non_normalised: Vector = (plane.point.0-ray.point.0, plane.point.1-ray.point.1, plane.point.2-ray.point.2);
    let t_not_normalised = normal_unit_vec.0*distances_non_normalised.0 + normal_unit_vec.1*distances_non_normalised.1 + normal_unit_vec.2*distances_non_normalised.2;
    let t: f64 = t_not_normalised/normalDotRayDir;

    if t > epsilon{
        return t;
    }
    return 0.0;
}


fn find_intersection_sphere(ray: Ray, sphere: Sphere) -> (f64, f64){
    let distance = (ray.point.0 - sphere.center.0, ray.point.1 - sphere.center.1, ray.point.2 - sphere.center.2);
    let a = ray.dir.0*ray.dir.0 + ray.dir.1*ray.dir.1 + ray.dir.2*ray.dir.2;
    let b = 2.0 * (distance.0*ray.dir.0 + distance.1*ray.dir.1 + distance.2*ray.dir.2);
    let c = (distance.0*distance.0+distance.1*distance.1+distance.2*distance.2) - (sphere.radius*sphere.radius);
    let delta = (b*b) - (4.0*a*c);

    if delta < 0.0{
        return (0.0, 0.0);
    }else if delta == 0.0{
        let point = - b / (2.0 * a);
        return (point, 0.0);
    }else{
        let t1 = (-b - delta.sqrt()) / (2.0 * a);
        let t2 = (-b + delta.sqrt()) / (2.0 * a);
        return (t1, t2);
    }
}


fn is_visible(t_to_check: f64, t_max: f64, ray: Ray, epsilon: f64, z_min: f64, z_max: f64) -> bool{
    let z_point = ray.point.2 + t_to_check * ray.dir.2;
    if (epsilon < t_to_check && t_to_check < t_max && z_min < z_point && z_point < z_max) {
        return true;
    }else {
        return false;
    }
}

fn find_closest_object(ray: Ray, objects: Vec<Plane>, epsilon:f64, z_min:f64, z_max:f64, spheres: Vec<Sphere>) -> (f64, Plane, f64, Sphere) {

    let c3: Point = (0.0,0.0,0.0);
    let rad3: f64 = 20.0;
    let col3: Color = (255.0,255.0,255.0);
    let k_diff3: f64 = 0.0;
    let k_refl3: f64 = 0.0;

    let falsesphere = Sphere {
        center:c3,
        radius:rad3,
        color: col3,
        k_diffuse:k_diff3,
        k_reflect:k_refl3,
    };


    let mut closestsphere = (999999999.9, falsesphere);

    for sphere in spheres{
        let sphere1 = sphere.clone();
        let ray1 = ray.clone();
        let t = find_intersection_sphere(ray1, sphere1);

        let t1 = t.0;
        let t2 = t.1;


        let ray1 = ray.clone();
        let ray2 = ray.clone();
        if t1 == 0.0{
            // nufink
        }else if is_visible(t1, closestsphere.0, ray1, epsilon, z_min, z_max) == true{
            closestsphere = (t1.clone(), sphere.clone());
        }else if is_visible(t2, closestsphere.0, ray2, epsilon, z_min, z_max) == true{
            closestsphere = (t2.clone(), sphere.clone());
        }
    }

    let p1: Point = (0.0,0.0,0.0);
    let n1: Normal = (0.0,0.0,0.0);
    let black: Color = (0.0, 0.0, 0.0);
    let diff1: f64 = 0.0;
    let refl1:f64 = 0.0;

    let falseplane = Plane {
        point:p1,
        normal:n1,
        color: black,
        diffuse:diff1,
        reflect:refl1,
    };

    let mut closest = (999999999.9, falseplane);

    for plane in objects{
        let plane1 = plane.clone();
        let ray1 = ray.clone();
        let t = find_intersection_plane(ray1, plane1, epsilon);

        let ray1 = ray.clone();
        if t == 0.0 {
            //nufink
        }else if is_visible(t, closest.0, ray1, epsilon, z_min, z_max) == true {
            let plane1 = plane.clone();
            closest = (t.clone(), plane1.clone());
        }
    }
    return (closest.0, closest.1, closestsphere.0, closestsphere.1);
}

fn reflect_ray(ray:Ray, surface_normal: Normal, point: Point, epsilon: f64) -> Ray{
    let twonormdotdir = 2.0*(ray.dir.0*surface_normal.0+ray.dir.1*surface_normal.1+ray.dir.2*surface_normal.2);
    let subtractor = (surface_normal.0*twonormdotdir, surface_normal.1*twonormdotdir, surface_normal.2*twonormdotdir);
    let new_dir = (ray.dir.0-subtractor.0, ray.dir.1-subtractor.1, ray.dir.2-subtractor.2);
    let new_dir_unit = convert_unit_dir(new_dir);

    let epsilonnorm = (epsilon*surface_normal.0, epsilon*surface_normal.1, epsilon*surface_normal.2);
    let pointepsmorm = (point.0+epsilonnorm.0, point.1+epsilonnorm.1, point.2+epsilonnorm.2);

    return Ray{
	point: pointepsmorm,
	dir: new_dir_unit,
    };
}

fn cast_ray(ray: Ray, depth: f64, objects:Vec<Plane>, spheres: Vec<Sphere>) -> Color {
    let epsilon = 0.003;
    let z_min:f64 = 200.0;
    let z_max:f64 = 1000.0;
    let light_position: Point = (500.0,500.0,500.0);
    let ambient_component = 0.05;
    let ambient_color: Color = (255.0,255.0,255.0);
    let max_depth: f64 = 2.0;


    let mut color: Color = (0.0,0.0,0.0);
    let ray1 = ray.clone();
    let objects1 = objects.clone();
    let spheres1 = spheres.clone();
    let intersection = find_closest_object(ray1, objects1, epsilon, z_min, z_max, spheres1);

    if intersection.0 < intersection.2 {
        if intersection.0 < 99999.9{
            let ray1 = ray.clone();
            let intersection_point = (ray1.point.0 + ray1.dir.0 * intersection.0, ray1.point.1 + ray1.dir.1 * intersection.0, ray1.point.2 + ray1.dir.2 * intersection.0);
            let light_dir: Dir = (light_position.0-intersection_point.0, light_position.1-intersection_point.1, light_position.2-intersection_point.2);
            let unit_light_dir = convert_unit_dir(light_dir);
            let light_ray = Ray{
                    point: intersection_point,
                    dir: unit_light_dir,
            };
            let objects1 = objects.clone();
            let spheres1 = spheres.clone();
            let light_intersection = find_closest_object(light_ray, objects1, epsilon, z_min, z_max, spheres1);
            let mut light: f64 = 1.0;
            if light_intersection.2 < light_intersection.0{ // fix later
                light = 0.0;
            }
            //colours etc
            let ambient_col: Color = (ambient_color.0*ambient_component*intersection.1.color.0/256.0, ambient_color.1*ambient_component*intersection.1.color.1/256.0, ambient_color.2*ambient_component*intersection.1.color.2/256.0);
            let unit_normal = convert_unit_normal(intersection.1.normal);
            let minval: f64 = 0.0;
            let lightandnorm = minval.max(unit_light_dir.0 * unit_normal.0 + unit_light_dir.1 * unit_normal.1 + unit_light_dir.2 * unit_normal.2);
            let diffuse_col: Color = (intersection.1.diffuse * lightandnorm * intersection.1.color.0, intersection.1.diffuse * lightandnorm * intersection.1.color.1, intersection.1.diffuse * lightandnorm * intersection.1.color.2);

            let new_ray = reflect_ray(ray, intersection.1.normal, intersection_point, epsilon);
            color = (ambient_col.0 + diffuse_col.0 * light, ambient_col.1 + diffuse_col.1 * light, ambient_col.2 + diffuse_col.2 * light);

            let objects1 = objects.clone();
            let spheres1 = spheres.clone();
            let new_ray_results = cast_ray(new_ray, depth + 1.0, objects1, spheres1);
            let new_ray_col: Color = (new_ray_results.0*intersection.1.reflect, new_ray_results.1*intersection.1.reflect, new_ray_results.2*intersection.1.reflect);
            let color_by_depth: Color = (new_ray_col.0/(depth+1.0), new_ray_col.1/(depth+1.0), new_ray_col.2/(depth+1.0));
            color = (color.0+color_by_depth.0, color.1+color_by_depth.1, color.2+color_by_depth.2);

        }
        return color;
    } else {
        if intersection.2 < 99999.9{
            let ray1 = ray.clone();
            let intersection_point = (ray1.point.0 + ray1.dir.0 * intersection.2, ray1.point.1 + ray1.dir.1 * intersection.2, ray1.point.2 + ray1.dir.2 * intersection.2);
            let light_dir: Dir = (light_position.0-intersection_point.0, light_position.1-intersection_point.1, light_position.2-intersection_point.2);
            let unit_light_dir = convert_unit_dir(light_dir);
            let light_ray = Ray{
                    point: intersection_point,
                    dir: unit_light_dir,
            };
            let objects1 = objects.clone();
            let spheres1 = spheres.clone();
            let light_intersection = find_closest_object(light_ray, objects1, epsilon, z_min, z_max, spheres1);
            let mut light: f64 = 1.0;
            if (light_intersection.3.center.0 != intersection.3.center.0) && light_intersection.2 < light_intersection.0{ // fix later
                light = 0.0;
            }
            //colours etc
            let ambient_col: Color = (ambient_color.0*ambient_component*intersection.3.color.0/256.0, ambient_color.1*ambient_component*intersection.3.color.1/256.0, ambient_color.2*ambient_component*intersection.3.color.2/256.0);
            let unit_normal = convert_unit_normal(intersection.3.getNormal(intersection_point));
            let minval: f64 = 0.0;
            let lightandnorm = minval.max(unit_light_dir.0 * unit_normal.0 + unit_light_dir.1 * unit_normal.1 + unit_light_dir.2 * unit_normal.2);
            let diffuse_col: Color = (intersection.3.k_diffuse * lightandnorm * intersection.3.color.0, intersection.3.k_diffuse * lightandnorm * intersection.3.color.1, intersection.3.k_diffuse * lightandnorm * intersection.3.color.2);

            let new_ray = reflect_ray(ray, intersection.3.getNormal(intersection_point), intersection_point, epsilon);
            color = (ambient_col.0 + diffuse_col.0 * light, ambient_col.1 + diffuse_col.1 * light, ambient_col.2 + diffuse_col.2 * light);

            let objects1 = objects.clone();
            let spheres1 = spheres.clone();
            let new_ray_results = cast_ray(new_ray, depth + 1.0, objects1, spheres1);
            let new_ray_col: Color = (new_ray_results.0*intersection.3.k_reflect, new_ray_results.1*intersection.3.k_reflect, new_ray_results.2*intersection.3.k_reflect);
            let color_by_depth: Color = (new_ray_col.0/(depth+1.0), new_ray_col.1/(depth+1.0), new_ray_col.2/(depth+1.0));
            color = (color.0+color_by_depth.0, color.1+color_by_depth.1, color.2+color_by_depth.2);

        }
        return color;
    }
}

fn generate_pixel_color(x:f64, y:f64, objects:Vec<Plane>, spheres: Vec<Sphere>) -> Color{
    let pixel_x = -49.95 + x * 0.1;
    let pixel_y = (999.0 - y) * 0.1 - 49.95;
    let pixel_z = 100.0;

    let dir: Dir = (pixel_x, pixel_y, pixel_z);
    let eye_point: Point = (0.0, 0.0, 0.0);
    let depth = 0.0;
    let initialRay = Ray{
        point: eye_point,
        dir: dir,
    };
    let spheres1 = spheres.clone();
    let color = cast_ray(initialRay, depth, objects, spheres1);

    return color;
}
fn generate_image(objects:Vec<Plane>, spheres: Vec<Sphere>) -> i32{
    let mut image = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(1000, 1000);
    for x in 0..=999 {
	println!("row {} is nomplete", x);
        for y in 0..=999 {
            let objects1 = objects.clone();
            let spheres1 = spheres.clone();
            let color = generate_pixel_color(x as f64, y as f64, objects1, spheres1);
            let r = cmp::min(color.0.round() as u8, 255);
            let g = cmp::min(color.1.round() as u8, 255);
            let b = cmp::min(color.2.round() as u8, 255);

            let pix = [r, g, b];

            image.put_pixel(x, y, Rgb(pix));
        }
    }

    image.save("output1.png");

    return 1;
}






fn main() {
    let mut objects: Vec<Plane> = Vec::new();
    let mut spheres: Vec<Sphere> = Vec::new();
    let epsilon = 0.003;
    let z_min = 200;
    let z_max = 1000;
    let light_position = vec![500,500,500];
    let ambient_component = 0.05;
    let ambient_color = vec![255,255,255];
    let max_depth = 2;


    let p1: Point = (0.0, -40.0, 0.0);
    let n1: Normal = (0.0,1.0,0.0);
    let red: Color = (120.0, 105.0, 186.0);
    let diff1: f64 = 0.5;
    let refl1:f64 = 1.0;

    let plane1 = Plane {
	point:p1,
	normal:n1,
	color: red,
	diffuse:diff1,
	reflect:refl1,
    };
    objects.push(plane1);

    let p2: Point = (0.0,0.0,999.0);
    let n2: Normal = (0.0,0.0,-1.0);
    let blue: Color = (10.0, 10.0, 10.0);
    let diff2: f64 = 1.0;
    let refl2:f64 = 0.5;

    let plane2 = Plane {
        point:p2,
        normal:n2,
        color: blue,
        diffuse:diff2,
        reflect:refl2,
    };
    objects.push(plane2);



    let c1: Point = (-30.0,-20.0,250.0);
    let rad1: f64 = 20.0;
    let col1: Color = (220.0,40.0,80.0);
    let k_diff1: f64 = (1.0);
    let k_refl1: f64 = (0.2);


    let sphere1 = Sphere {
        center:c1,
        radius:rad1,
        color: col1,
        k_diffuse:k_diff1,
        k_reflect:k_refl1,
    };
    spheres.push(sphere1);

    let c2: Point = (-30.0,-20.0,250.0);
    let rad2: f64 = 20.0;
    let col2: Color = (70.0,30.0,255.0);
    let k_diff2: f64 = (1.0);
    let k_refl2: f64 = (0.05);

    let sphere2 = Sphere {
        center:c2,
        radius:rad2,
        color: col2,
        k_diffuse:k_diff2,
        k_reflect:k_refl2,
    };
    spheres.push(sphere2);

    println!("inputs done");

    generate_image(objects, spheres);

}

