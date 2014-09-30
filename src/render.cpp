// Source file for photonmap renderer
// This implementation is a simple raycaster.
// Replace it with your own code.



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Graphics/R3Graphics.h"


////////////////////////////////////////////////////////////////////////
// Constant power photon of a specific color
////////////////////////////////////////////////////////////////////////

class Photon {
public:
  enum { R, G, B };

  Photon(R3Ray ray, int color);
  void Draw(double radius) const;

  const R3Ray ray;
  const R3Vector from_dir;
  const int color;
};

Photon::Photon(R3Ray ray, int color)
  : ray(ray), color(color)
{
}

void
Photon::Draw(double radius) const
{
  R3Point start = this->ray.Start();
  R3Vector dir = this->ray.Vector();

  if (color == R) glColor3d(1, 0, 0);
  if (color == G) glColor3d(0, 1, 0);
  if (color == B) glColor3d(0, 0, 1);

  R3Span(start, start + dir * radius).Draw();
}

////////////////////////////////////////////////////////////////////////
// Randomly sample light sources based on their intensity.
// Generate photons from light sources.
////////////////////////////////////////////////////////////////////////

RNScalar Random() {
  return static_cast <RNScalar> (rand())
       / static_cast <RNScalar> (RAND_MAX);
}

R3Vector RandomVectorUniform() {
  while (true) {
    RNScalar x = 2.0 * Random() - 1.0;
    RNScalar y = 2.0 * Random() - 1.0;
    RNScalar z = 2.0 * Random() - 1.0;

    if (x*x + y*y + z*z > 1.0) continue;

    R3Vector v = R3Vector(x, y, z);
    v.Normalize();
    return v;
  }
}

R3Vector RandomVectorInDir(R3Vector dir) {
  while (true) {
    R3Vector rand = RandomVectorUniform();
    if (rand.Dot(dir) > 0) return rand;
  }
}

int PickColor(R3Light *light) {
  RNRgb color = light->Color();
  RNScalar total = color.R() + color.G() + color.B();
  RNScalar r = Random() * total;

  if (r < color.R()) return Photon::R;
  r -= color.R();
  if (r < color.G()) return Photon::G;
  r -= color.G();
  if (r < color.B()) return Photon::B;

  printf("ERROR: color wat\n");
  return -1;
}

R3Vector RandomVectorSpot(R3SpotLight *light) {
  RNAngle cutoffangle = light->CutOffAngle();
  RNScalar dropoffrate = light->DropOffRate();
  R3Vector dir = light->Direction();

  while (true) {
    R3Vector v = RandomVectorUniform();
    RNScalar cos_alpha = v.Dot(dir);
    if (cos_alpha < 0) continue;
    if (cos(cutoffangle) > cos_alpha) continue;

    RNScalar r = Random();
    if (r < pow(cos_alpha, dropoffrate)) return v;
  }
}

Photon *
PhotonFromDirLight(R3DirectionalLight *light, int scene_radius)
{
  R3Vector dir = light->Direction();
  dir.Normalize();

  R3Vector norm = R3Vector(0,1,0);

  R3Vector axis = R3Vector(norm);
  axis.Cross(dir);
  axis.Normalize();

  RNScalar cos_angle = dir.Dot(norm);
  RNAngle angle = acos(cos_angle);

  RNScalar x, z;
  while (true) {
    x = 1.5 * scene_radius * (2.0 * Random() - 1.0);
    z = 1.5 * scene_radius * (2.0 * Random() - 1.0);
    if (x*x + z*z < 2 * scene_radius * scene_radius) break;
  }
  R3Vector pos = R3Vector(x,0,z);

  pos.Rotate(axis, angle);
  pos -= scene_radius * 2 * dir;

  R3Ray ray = R3Ray(pos.Point(), dir);
  Photon *photon = new Photon(ray, PickColor(light));

  return photon;
}
Photon *
PhotonFromPointLight(R3PointLight *light)
{
  R3Ray ray = R3Ray(light->Position(), RandomVectorUniform());
  Photon *photon = new Photon(ray, PickColor(light));

  return photon;
}
Photon *
PhotonFromSpotLight(R3SpotLight *light)
{
  R3Ray ray = R3Ray(light->Position(), RandomVectorSpot(light));
  Photon *photon = new Photon(ray, PickColor(light));

  return photon;
}


RNScalar TotalLightIntensity(R3Scene *scene) {
  RNScalar total = 0;
  for (int i = 0; i < scene->NLights(); i ++) {
    R3Light *light = scene->Light(i);
    total += light->Intensity();
  }
  return total;
}

R3Light *RandomLight(R3Scene *scene, RNScalar total_intensity) {
  RNScalar r = total_intensity * Random();

  for (int i = 0; i < scene->NLights(); i ++) {
    R3Light *light = scene->Light(i);
    r -= light->Intensity();

    if (r < 0) {
      return light;
    }
  }
  return NULL;
}

RNArray<Photon *> *
PhotonsFromLights(R3Scene *scene, int num)
{
  printf("Generating photons from lights.\n");
  RNArray<Photon *> *photons = new RNArray<Photon *>;

  RNScalar total_intensity = TotalLightIntensity(scene);
  printf("Total light intensity %f.\n", total_intensity);

  double radius = scene->BBox().DiagonalRadius();
  for (int i = 0; i < num; i ++) {
    R3Light *light = RandomLight(scene, total_intensity);

    int light_class = light->ClassID();

    if (light_class == R3DirectionalLight::CLASS_ID()) {
      photons->Insert(PhotonFromDirLight((R3DirectionalLight *) light, radius));
    }
    else if (light_class == R3PointLight::CLASS_ID()) {
      photons->Insert(PhotonFromPointLight((R3PointLight *) light));
    }
    else if (light_class == R3SpotLight::CLASS_ID()) {
      photons->Insert(PhotonFromSpotLight((R3SpotLight *) light));
    }
  }

  return photons;
}

////////////////////////////////////////////////////////////////////////
// Function for reflections and transmissions
////////////////////////////////////////////////////////////////////////

Photon *
DiffuseBounce(Photon *source_photon, R3Point pos, R3Vector norm) {
  R3Vector dir = RandomVectorInDir(norm);
  R3Ray ray = R3Ray(pos, dir);
  Photon *photon = new Photon(ray, source_photon->color);
  return photon;
}

Photon *
SpecularBounce(Photon *source_photon,
    R3Point inters_pos, R3Vector inters_norm) {
  R3Vector source_dir = -source_photon->ray.Vector();
  R3Vector dir = 2 * inters_norm * inters_norm.Dot(source_dir) - source_dir;

  R3Ray ray = R3Ray(inters_pos, dir);
  Photon *photon = new Photon(ray, source_photon->color);
  return photon;
}

Photon *
TransmissionBounce(Photon *source_photon,
    RNScalar index_of_refraction,
    R3Point inters_pos, R3Vector inters_norm) {
  R3Vector source_dir = -source_photon->ray.Vector();
  inters_norm.Normalize();
  source_dir.Normalize();

  RNScalar ir_ratio;
  if (source_dir.Dot(inters_norm) > 0) {
    ir_ratio = 1.0 / index_of_refraction;
  } else {
    ir_ratio = index_of_refraction / 1.0;
  }

  // Snell's law
  RNScalar theta_i = acos(source_dir.Dot(inters_norm));
  RNScalar theta_r = asin(ir_ratio * sin(theta_i));

  R3Vector dir = inters_norm * (ir_ratio * cos(theta_i) - cos(theta_r))
                 - source_dir * ir_ratio;

  R3Ray ray = R3Ray(inters_pos, dir);
  Photon *photon = new Photon(ray, source_photon->color);
  return photon;
}

Photon *
ScatterPhoton(Photon *source_photon, R3Scene *scene) {
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  RNBoolean intersected = scene->Intersects(source_photon->ray,
    &node, &element, &shape, &point, &normal, &t);

  if (intersected) {
    R3Material *material = element->Material();
    const R3Brdf *brdf = material->Brdf();
    int color = source_photon->color;

    RNScalar kd = brdf->Diffuse()[color];
    RNScalar ks = brdf->Specular()[color];
    RNScalar kt = brdf->Transmission()[color];

    RNScalar r = Random();
    if (r < kd) {
      // diffuse bounce
      return DiffuseBounce(source_photon, point, normal);
    } else if (r < ks + kd) {
      // specular bounce
      return SpecularBounce(source_photon, point, normal);
    } else if (r < kt + ks + kd) {
      // transmission
      RNScalar ir = brdf->IndexOfRefraction();
      return TransmissionBounce(source_photon, ir, point, normal);
    }
    r -= kt;
  }

  return NULL;
}

RNArray<Photon *> *
ScatterPhotons(RNArray<Photon *> *source_photons, R3Scene *scene)
{
  RNArray<Photon *> *scattered_photons = new RNArray<Photon *>;

  for (int i = 0; i < source_photons->NEntries(); i ++) {
    Photon *source_photon = source_photons->Kth(i);
    Photon *scattered_photon = ScatterPhoton(source_photon, scene);
    if (scattered_photon != NULL) {
      scattered_photons->Insert(scattered_photon);
    }
  }

  return scattered_photons;
}



////////////////////////////////////////////////////////////////////////
// Function to draw photons for debugging
////////////////////////////////////////////////////////////////////////

static RNArray<Photon *> *cached_photons = NULL;

RNArray<Photon *> *
GenerateAndCachePhotons(R3Scene *scene) {
  if (cached_photons != NULL) {
    return cached_photons;
  }
  RNArray<Photon *> *photons1 = PhotonsFromLights(scene, 100000);
  RNArray<Photon *> *photons2 = ScatterPhotons(photons1, scene);
  RNArray<Photon *> *photons3 = ScatterPhotons(photons2, scene);

  cached_photons = new RNArray<Photon *>;
  cached_photons->Append(*photons1);
  cached_photons->Append(*photons2);
  cached_photons->Append(*photons3);

  return cached_photons;
}

void
DrawPhotons(R3Scene *scene)
{
  double radius = 0.025 * scene->BBox().DiagonalRadius();
  RNArray<Photon *> *photons = GenerateAndCachePhotons(scene);
  for (int i = 0; i < photons->NEntries(); i ++) {
    Photon *photon = photons->Kth(i);
    if (photon != NULL) {
      photon->Draw(radius);
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Function to render image with photon mapping
////////////////////////////////////////////////////////////////////////

R2Image *
RenderImage(R3Scene *scene,
  int width, int height,
  int print_verbose)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int ray_count = 0;

  // Allocate image
  R2Image *image = new R2Image(width, height);
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    return NULL;
  }

  // Convenient variables
  const R3Point& eye = scene->Camera().Origin();
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  // Draw intersection point and normal for some rays
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      R3Ray ray = scene->Viewer().WorldRay(i, j);
      if (scene->Intersects(ray, &node, &element, &shape, &point, &normal, &t)) {
        // Get intersection information
        const R3Material *material = (element) ? element->Material() : NULL;
        const R3Brdf *brdf = (material) ? material->Brdf() : NULL;

        // Compute color
        RNRgb color = scene->Ambient();
        if (brdf) {
          color += brdf->Emission();
          for (int k = 0; k < scene->NLights(); k++) {
            R3Light *light = scene->Light(k);
            color += light->Reflection(*brdf, eye, point, normal);
          }
        }

        // Set pixel color
        image->SetPixelRGB(i, j, color);

        // Update ray count
        ray_count++;
      }
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Rendered image ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Rays = %d\n", ray_count);
    fflush(stdout);
  }

  // Return image
  return image;
}
