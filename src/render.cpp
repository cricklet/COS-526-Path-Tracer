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
// Photon samples. These provide radiant flux information.
////////////////////////////////////////////////////////////////////////

struct PhotonSample {
  R3Point position;
  R3Vector incident;
  RNRgb rgb;
};

static R3Point
GetPhotonSamplePosition(PhotonSample *photon, void *dummy)
{
  // Return point position (used for Kdtree)
  return photon->position;
}

void DrawPhotonSample(PhotonSample *photon, double radius) {
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
  enum { LIGHT, DIFFUSE, SPECULAR, TRANSMISSION };
  enum { GLOBAL, CAUSTIC };

  Photon(R3Ray ray, int color, int photon_type);
  Photon(R3Ray ray, Photon *photon, int photon_type);
  void Draw(double radius) const;
  void DrawPath(double radius) const;

  const R3Ray ray;
  const Photon *source;
  const int color;
  const int photon_type;
  const int path_type;

private:
  int GetPathType();
};

int
Photon::GetPathType()
{
  if (photon_type == DIFFUSE && source != NULL) {
    int source_type = source->photon_type;
    if (source_type == SPECULAR || source_type == TRANSMISSION ) {
      return CAUSTIC;
    }
  }

  return GLOBAL;
}

Photon::Photon(R3Ray ray, Photon *source, int photon_type)
  : ray(ray), source(source), color(source->color),
    photon_type(photon_type), path_type(this->GetPathType())
{
}

Photon::Photon(R3Ray ray, int color, int photon_type)
  : ray(ray), source(NULL), color(color),
    photon_type(photon_type), path_type(this->GetPathType())
{
}

void
Photon::Draw(double radius) const
{
  if (color == R) glColor3d(1,0,0);
  if (color == G) glColor3d(0,1,0);
  if (color == B) glColor3d(0,0,1);

  R3Point start = this->ray.Start();
  R3Vector dir = this->ray.Vector();

  R3Span(start, start + dir * radius).Draw();
}

void
Photon::DrawPath(double radius) const
{
  if (color == R) glColor3d(1,0,0);
  if (color == G) glColor3d(0,1,0);
  if (color == B) glColor3d(0,0,1);

  const Photon *photon = this;
  while (photon != NULL && photon->source != NULL) {
    R3Point start = photon->ray.Start();
    R3Point end = photon->source->ray.Start();

    R3Span(start, end).Draw();

    photon = photon->source;
  }
}

////////////////////////////////////////////////////////////////////////
// Randomly sample light sources based on their intensity.
// Generate photons from light sources.
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
  Photon *photon = new Photon(ray, PickColor(light), Photon::LIGHT);

  return photon;
}
Photon *
PhotonFromPointLight(R3PointLight *light)
{
  R3Ray ray = R3Ray(light->Position(), RandomVectorUniform());
  Photon *photon = new Photon(ray, PickColor(light), Photon::LIGHT);

  return photon;
}
Photon *
PhotonFromSpotLight(R3SpotLight *light)
{
  R3Ray ray = R3Ray(light->Position(), RandomVectorSpot(light));
  Photon *photon = new Photon(ray, PickColor(light), Photon::LIGHT);

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
  Photon *photon = new Photon(ray, source_photon, Photon::DIFFUSE);
  return photon;
}

Photon *
SpecularBounce(Photon *source_photon,
    R3Point inters_pos, R3Vector inters_norm) {
  R3Vector source_dir = -source_photon->ray.Vector();
  R3Vector dir = 2 * inters_norm * inters_norm.Dot(source_dir) - source_dir;

  R3Ray ray = R3Ray(inters_pos, dir);
  Photon *photon = new Photon(ray, source_photon, Photon::SPECULAR);
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
  Photon *photon = new Photon(ray, source_photon, Photon::TRANSMISSION);
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

PhotonSample *
ProcessPhotonSample(Photon *photon) {
  if (photon->source == NULL) return NULL;

  R3Point position = photon->ray.Start();
  R3Point source = photon->source->ray.Start();
  R3Vector incident = position - source;
  incident.Normalize();

  PhotonSample *photon_sample = new PhotonSample;
  photon_sample->position = position;
  photon_sample->incident = incident;

  int color = photon->color;
  RNRgb rgb;
  if (color == Photon::R)      rgb = RNRgb(1,0,0);
  else if (color == Photon::G) rgb = RNRgb(0,1,0);
  else if (color == Photon::B) rgb = RNRgb(0,0,1);

  photon_sample->rgb = rgb;

  return photon_sample;
}

RNArray<PhotonSample *> *
GeneratePhotonSamples(RNArray<Photon *> *photons, const R3Box& bbox) {
  RNArray<PhotonSample *> *photon_samples = new RNArray<PhotonSample *>;

  for (int i = 0; i < photons->NEntries(); i ++) {
    Photon *photon = photons->Kth(i);
    PhotonSample *photon_sample = ProcessPhotonSample(photon);
    if (photon_sample != NULL) {
      photon_samples->Insert(photon_sample);
    }
  }

  return photon_samples;
}

////////////////////////////////////////////////////////////////////////
// Function for generating + caching photons
// (both path tracing + rendering)
////////////////////////////////////////////////////////////////////////

static RNArray<Photon *> *cached_all_photons = NULL;
static RNArray<Photon *> *cached_final_photons = NULL;

static RNArray<PhotonSample *> *cached_global_samples = NULL;
static R3Kdtree<PhotonSample *> *cached_kd_global_samples = NULL;

static RNArray<PhotonSample *> *cached_caustic_samples = NULL;
static R3Kdtree<PhotonSample *> *cached_kd_caustic_samples = NULL;

void
CachePhotons(R3Scene *scene) {
  if (cached_final_photons  != NULL
      && cached_all_photons != NULL) {
    return;
  }

  int initial_photons = 100000;
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

  cached_all_photons = photons;
  cached_final_photons = final_photons;
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

void
CacheSamples(R3Scene *scene) {
  if (cached_global_samples     != NULL
      && cached_caustic_samples != NULL) {
    return;
  }
  printf("Generating photon samples\n");

  RNArray<Photon *> *photons = GetAllPhotons(scene);

  RNArray<Photon *> *caustic_photons = new RNArray<Photon *>;
  RNArray<Photon *> *global_photons = new RNArray<Photon *>;

  for (int i = 0; i < photons->NEntries(); i ++) {
    Photon *photon = photons->Kth(i);
    int type = photon->path_type;
    if (type == Photon::GLOBAL) {
      global_photons->Insert(photon);
    } else if (type == Photon::CAUSTIC){
      caustic_photons->Insert(photon);
    }
  }

  RNArray<PhotonSample *> *global_samples =
    GeneratePhotonSamples(global_photons, scene->BBox());
  RNArray<PhotonSample *> *caustic_samples =
    GeneratePhotonSamples(caustic_photons, scene->BBox());

  R3Kdtree<PhotonSample *> *kd_global_samples =
    new R3Kdtree<PhotonSample *>(*global_samples, GetPhotonSamplePosition);
  R3Kdtree<PhotonSample *> *kd_caustic_samples =
    new R3Kdtree<PhotonSample *>(*caustic_samples, GetPhotonSamplePosition);

  delete caustic_photons;
  delete global_photons;

  printf("Stored photon samples into kd trees\n");

  cached_global_samples = global_samples;
  cached_kd_global_samples = kd_global_samples;

  cached_caustic_samples = caustic_samples;
  cached_kd_caustic_samples = kd_caustic_samples;
}

RNArray<PhotonSample *> *
GetGlobalSamples(R3Scene *scene) {
  CacheSamples(scene);
  return cached_global_samples;
}

RNArray<PhotonSample *> *
GetCausticSamples(R3Scene *scene) {
  CacheSamples(scene);
  return cached_caustic_samples;
}

R3Kdtree<PhotonSample *> *
GetKdGlobalSamples(R3Scene *scene) {
  CacheSamples(scene);
  return cached_kd_global_samples;
}

R3Kdtree<PhotonSample *> *
GetKdCausticSamples(R3Scene *scene) {
  CacheSamples(scene);
  return cached_kd_caustic_samples;
}

////////////////////////////////////////////////////////////////////////
// Function to draw debugging data
////////////////////////////////////////////////////////////////////////

void
DrawCausticSamples(R3Scene *scene)
{
  double radius = 0.01 * scene->BBox().DiagonalRadius();

  RNArray<PhotonSample *> *photons = GetCausticSamples(scene);
  for (int i = 0; i < photons->NEntries(); i ++) {
    DrawPhotonSample(photons->Kth(i), radius);
  }
}

void
DrawGlobalSamples(R3Scene *scene)
{
  double radius = scene->BBox().DiagonalRadius();

  RNArray<PhotonSample *> *photons = GetGlobalSamples(scene);
  for (int i = 0; i < photons->NEntries(); i ++) {
    DrawPhotonSample(photons->Kth(i), 0.01 * radius);
  }
}

void
DrawPhotons(R3Scene *scene)
{
  double radius = scene->BBox().DiagonalRadius();

  RNArray<Photon *> *all = GetAllPhotons(scene);
  for (int i = 0; i < all->NEntries(); i ++) {
    Photon *photon = all->Kth(i);
    photon->Draw(0.01 * radius);
  }

  RNArray<Photon *> *final = GetAllPhotons(scene);
  for (int i = 0; i < final->NEntries(); i ++) {
    Photon *photon = all->Kth(i);
    if (Random() > 0.999) {
      photon->DrawPath(0.01 * radius);
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Function to render image with photon mapping
////////////////////////////////////////////////////////////////////////

RNRgb RenderPixel(R3Scene *scene, R3Ray ray) {
  double radius = scene->BBox().DiagonalRadius();

  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  RNBoolean intersection = scene->Intersects(ray,
    &node, &element, &shape, &point, &normal, &t);

  if (!intersection) {
    return scene->Ambient();
  }

  R3Kdtree<PhotonSample *> *global_samples = GetKdGlobalSamples(scene);
  R3Kdtree<PhotonSample *> *caustic_samples = GetKdCausticSamples(scene);

  RNArray<PhotonSample *> samples = RNArray<PhotonSample *>();

  global_samples->FindClosest(point, 0, radius * 0.05, 500, samples);
  caustic_samples->FindClosest(point, 0, radius * 0.01, 500, samples);

  RNRgb color = RNRgb(0,0,0);

  // If we assume each of these flux samples is passing through the
  // intersection point, how much of that flux is reflected towards the viewer?
  // R3Material *material = element->Material();
  // const R3Brdf *brdf = material->Brdf();
  // RNScalar kd = brdf->Diffuse()[color];
  // RNScalar ks = brdf->Specular()[color];
  // RNScalar kt = brdf->Transmission()[color];
  // RNScalar k_total = fmax(1.0, kd + ks + kt);
  //
  // R3Vector viewer_dir = ray.Start() - point;
  //
  // RNScalar num_samples = samples.NEntries();
  // for (int i = 0; i < num_samples; i ++) {
  //   PhotonSample *sample = samples.Kth(i);
  //
  //   R3Vector incident_dir = sample->incident;
  // }

  return color;
}

R2Image *
RenderImage(R3Scene *scene,
  int width, int height,
  int print_verbose)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate image
  R2Image *image = new R2Image(width, height);
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    return NULL;
  }

  for (int i = 0; i < width; i++) {
    printf("Rendering pixel col %d out of %d\n", i, width);
    for (int j = 0; j < height; j++) {
      R3Ray ray = scene->Viewer().WorldRay(i, j);
      RNRgb color = RenderPixel(scene, ray);
      image->SetPixelRGB(i, j, color);
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Rendered image ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return image
  return image;
}
