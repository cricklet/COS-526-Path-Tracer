// Source file for photonmap renderer
// This implementation is a simple raycaster.
// Replace it with your own code.



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Graphics/R3Graphics.h"


RNScalar Random() {
  return static_cast <RNScalar> (rand())
       / static_cast <RNScalar> (RAND_MAX);
}

////////////////////////////////////////////////////////////////////////
// Constant power photon of a specific color for path tracing
////////////////////////////////////////////////////////////////////////

struct RenderPhoton {
  R3Point position;
  R3Vector incident;
  RNRgb rgb;
};

static R3Point
GetRenderPhotonPointPosition(RenderPhoton *photon, void *dummy)
{
  // Return point position (used for Kdtree)
  return photon->position;
}

void DrawRenderPhoton(RenderPhoton *photon, double radius) {
  glColor3d(photon->rgb.R(), photon->rgb.G(), photon->rgb.B());

  R3Point pos = photon->position;
  R3Span(pos, pos - photon->incident * radius).Draw();
}

////////////////////////////////////////////////////////////////////////
// Constant power photon of a specific color for path tracing
////////////////////////////////////////////////////////////////////////

class Photon {
public:
  enum { R, G, B };

  Photon(R3Ray ray, int color);
  Photon(R3Ray ray, Photon *photon);
  void Draw(double radius) const;
  void DrawPath(double radius) const;

  R3Point position;
  R3Vector incident;

  R3Vector normal;
  R3Material *material;

  R3Point source;
  const int color;

  const Photon *parent;
};

Photon::Photon(R3Point position, Photon *parent)
  : position(position),
    incident(position - parent->position),
    source(parent->position),
    color(parent->color),
    parent(parent)
{
  incident.Normalize();
}

Photon::Photon(R3Point position, R3Point source, int color)
  : position(position),
    incident(position - source),
    source(source),
    color(color),
    parent(NULL)
{
  incident.Normalize();
}

void
Photon::Draw(double radius) const
{
  if (color == R) glColor3d(1,0,0);
  if (color == G) glColor3d(0,1,0);
  if (color == B) glColor3d(0,0,1);

  R3Span(position, position + incident * radius).Draw();
}

void
Photon::DrawPath(double radius) const
{
  if (color == R) glColor3d(1,0,0);
  if (color == G) glColor3d(0,1,0);
  if (color == B) glColor3d(0,0,1);

  const Photon *photon = this;
  while (photon != NULL) {
    R3Span(photon->position, photon->position).Draw();
    photon = photon->parent;
  }
}

////////////////////////////////////////////////////////////////////////
// Generate rays from lights to shoot photons.
////////////////////////////////////////////////////////////////////////

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

R3Ray
RayFromDirLight(R3DirectionalLight *light, int scene_radius)
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

  return R3Ray(pos.Point(), dir);
}

R3Ray
RayFromPointLight(R3PointLight *light)
{
  return R3Ray(light->Position(), RandomVectorUniform());
}

R3Ray
RayFromSpotLight(R3SpotLight *light)
{
  return R3Ray(light->Position(), RandomVectorSpot(light));
}

////////////////////////////////////////////////////////////////////////
// Generate photons.
////////////////////////////////////////////////////////////////////////

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

Photon
ShootPhoton(R3Ray &ray, int color, R3Scene *scene, Photon *parent=NULL) {
  return Photon(
    ray.Start(),
  )
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
  Photon *photon = new Photon(ray, source_photon);
  return photon;
}

Photon *
SpecularBounce(Photon *source_photon,
    R3Point inters_pos, R3Vector inters_norm) {
  R3Vector source_dir = -source_photon->ray.Vector();
  R3Vector dir = 2 * inters_norm * inters_norm.Dot(source_dir) - source_dir;

  R3Ray ray = R3Ray(inters_pos, dir);
  Photon *photon = new Photon(ray, source_photon);
  return photon;
}

Photon *
TransmissionBounce(Photon *source_photon,
    RNScalar index_of_refraction,
    R3Point inters_pos, R3Vector inters_norm) {

  R3Vector l = - source_photon->ray.Vector();
  R3Vector n;
  RNScalar ni;
  RNScalar nr;

  if (l.Dot(inters_norm) > 0) {
    ni = 1.0;
    nr = index_of_refraction;
    n = inters_norm;
  } else {
    ni = index_of_refraction;
    nr = 1.0;
    n = - inters_norm;
  }

  RNScalar ratio = ni / nr;

  // Snell's law
  RNScalar theta_i = acos(l.Dot(n));
  RNScalar theta_r = asin(ratio * sin(theta_i));

  R3Vector dir = n * (ratio * cos(theta_i) - cos(theta_r))
               - l * ratio;

  R3Ray ray = R3Ray(inters_pos, dir);
  Photon *photon = new Photon(ray, source_photon);
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

  // offset the ray a bit
  R3Ray ray = R3Ray(source_photon->ray.Point(0.01),
                    source_photon->ray.Vector());

  RNBoolean intersected = scene->Intersects(ray,
    &node, &element, &shape, &point, &normal, &t);

  if (intersected) {
    R3Material *material = element->Material();
    const R3Brdf *brdf = material->Brdf();
    int color = source_photon->color;

    RNScalar kd = brdf->Diffuse()[color];
    RNScalar ks = brdf->Specular()[color];
    RNScalar kt = brdf->Transmission()[color];
    RNScalar k_total = fmax(1.0, kd + ks + kt);

    RNScalar r = Random() * k_total;
    if (r < kd) {
      return DiffuseBounce(source_photon, point, normal);
    } else if (r < ks + kd) {
      return SpecularBounce(source_photon, point, normal);
    } else if (r < kt + ks + kd) {
      RNScalar ir = brdf->IndexOfRefraction();
      return TransmissionBounce(source_photon, ir, point, normal);
    }
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
// Store the photons in a KDTree and render them!
////////////////////////////////////////////////////////////////////////

RenderPhoton *
ProcessRenderPhoton(Photon *photon) {
  if (photon->source == NULL) return NULL;

  R3Point position = photon->ray.Start();
  R3Point source = photon->source->ray.Start();
  R3Vector incident = position - source;
  incident.Normalize();

  RenderPhoton *render_photon = new RenderPhoton;
  render_photon->position = position;
  render_photon->incident = incident;

  int color = photon->color;
  RNRgb rgb;
  if (color == Photon::R)      rgb = RNRgb(1,0,0);
  else if (color == Photon::G) rgb = RNRgb(0,1,0);
  else if (color == Photon::B) rgb = RNRgb(0,0,1);

  render_photon->rgb = rgb;

  return render_photon;
}

RNArray<RenderPhoton *> *
GenerateRenderPhotons(RNArray<Photon *> *photons, const R3Box& bbox) {
  RNArray<RenderPhoton *> *render_photons = new RNArray<RenderPhoton *>;

  for (int i = 0; i < photons->NEntries(); i ++) {
    Photon *photon = photons->Kth(i);
    RenderPhoton *render_photon = ProcessRenderPhoton(photon);
    if (render_photon != NULL) {
      render_photons->Insert(render_photon);
    }
  }

  return render_photons;
}

////////////////////////////////////////////////////////////////////////
// Function for generating + caching photons
// (both path tracing + rendering)
////////////////////////////////////////////////////////////////////////

static RNArray<Photon *> *cached_all_photons = NULL;
static RNArray<Photon *> *cached_final_photons = NULL;

static RNArray<RenderPhoton *> *cached_render_photons = NULL;
static R3Kdtree<RenderPhoton *> *cached_kd_render_photons = NULL;

void
CachePhotons(R3Scene *scene) {
  if (cached_final_photons     != NULL
      && cached_render_photons != NULL
      && cached_all_photons    != NULL) {
    return;
  }

  int initial_photons = 10000;
  int total_photons = initial_photons * 6;

  RNArray<Photon *> *photons = PhotonsFromLights(scene, initial_photons);
  printf("Creating initial %d photons\n", photons->NEntries());

  RNArray<Photon *> *final_photons = new RNArray<Photon *>;
  for (int i = 0; i < photons->NEntries(); i ++) {
    Photon *photon = photons->Kth(i);
    Photon *scattered_photon = ScatterPhoton(photon, scene);
    if (scattered_photon != NULL) {
      photons->Insert(scattered_photon);
    } else {
      final_photons->Insert(photon);
    }

    if (photons->NEntries() > total_photons) {
      break;
    }
  }

  printf("Finished with %d photons\n", photons->NEntries());

  RNArray<RenderPhoton *> *render_photons =
    GenerateRenderPhotons(photons, scene->BBox());
  R3Kdtree<RenderPhoton *> *kd_render_photons =
    new R3Kdtree<RenderPhoton *>(*render_photons, GetRenderPhotonPointPosition);

  printf("Moved %d render photons into a kd tree\n",
         render_photons->NEntries());

  cached_all_photons = photons;
  cached_final_photons = final_photons;
  cached_render_photons = render_photons;
  cached_kd_render_photons = kd_render_photons;
}

RNArray<Photon *> *
GetFinalPhotons(R3Scene *scene) {
  CachePhotons(scene);
  return cached_final_photons;
}

RNArray<Photon *> *
GetAllPhotons(R3Scene *scene) {
  CachePhotons(scene);
  return cached_all_photons;
}

RNArray<RenderPhoton *> *
GetRenderPhotons(R3Scene *scene) {
  CachePhotons(scene);
  return cached_render_photons;
}

////////////////////////////////////////////////////////////////////////
// Function to draw debugging data
////////////////////////////////////////////////////////////////////////

void
DrawRenderPhotons(R3Scene *scene)
{
  double radius = 0.025 * scene->BBox().DiagonalRadius();

  RNArray<RenderPhoton *> *photons = GetRenderPhotons(scene);
  for (int i = 0; i < photons->NEntries(); i ++) {
    DrawRenderPhoton(photons->Kth(i), radius);
  }
}

void
DrawPhotons(R3Scene *scene)
{
  double radius = 0.01 * scene->BBox().DiagonalRadius();

  RNArray<Photon *> *all = GetAllPhotons(scene);
  for (int i = 0; i < all->NEntries(); i ++) {
    Photon *photon = all->Kth(i);
    photon->Draw(radius);
  }

  RNArray<Photon *> *final = GetAllPhotons(scene);
  for (int i = 0; i < final->NEntries(); i ++) {
    Photon *photon = all->Kth(i);
    if (Random() > 0.999) {
      photon->DrawPath(radius);
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
