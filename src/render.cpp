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

  glColor3d(1.0, 1.0, 1.0);
  R3Span(start, start + dir * 2 * radius).Draw();
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

Photon *
PhotonFromDirLight(R3DirectionalLight *light, int scene_radius)
{
  return NULL;
}
Photon *
PhotonFromPointLight(R3PointLight *light)
{
  R3Ray ray = R3Ray(light->Position(), RandomVectorUniform());
  Photon *photon = new Photon(ray, Photon::R);

  return photon;
}
Photon *
PhotonFromSpotLight(R3SpotLight *light)
{
  return NULL;
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

static RNArray<Photon *> *cached_light_photons = NULL;

RNArray<Photon *> *
PhotonsFromLights(R3Scene *scene, int num)
{
  // Memoize the result
  if (cached_light_photons != NULL) {
    return cached_light_photons;
  }

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

  cached_light_photons = photons;
  return photons;
}

////////////////////////////////////////////////////////////////////////
// Function to draw photons for debugging
////////////////////////////////////////////////////////////////////////

void
DrawPhotons(R3Scene *scene)
{
  double radius = 0.025 * scene->BBox().DiagonalRadius();

  // Draw photons coming out of light sources
  RNArray<Photon *> *photons = PhotonsFromLights(scene, 100);
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
