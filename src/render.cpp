// Source file for photonmap renderer
// This implementation is a simple raycaster.
// Replace it with your own code.



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Graphics/R3Graphics.h"


RNRgb
ClampRGB(RNRgb c) {
  return RNRgb(fmin(1.0, c.R()), fmin(1.0, c.G()), fmin(1.0, c.B()));
}

RNRgb
ScaleRGB(RNRgb c) {
  RNScalar c_max = fmax(fmax(c.R(), c.G()), c.B());
  if (c_max < 1.0) return c;
  return c / c_max;
}

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
// Constant power photon for path tracing
////////////////////////////////////////////////////////////////////////

class Path {
public:
  enum { LIGHT, DIFFUSE, SPECULAR, TRANSMISSION };

  Photon(R3Ray ray, RNRgb rgb, int photon_type, Photon *source);

  void Draw(double radius) const;
  void DrawPath(double radius) const;

  const R3Point start;
  const R3Point end;
  const Path *previous;
  const RNRgb rgb;
  const int start_type;
};

int
GetPathType(int photon_type, Photon *source)
{
  if (photon_type == Photon::DIFFUSE && source != NULL) {
    int source_type = source->photon_type;
    if (source_type == Photon::SPECULAR
        || source_type == Photon::TRANSMISSION ) {
      return Photon::CAUSTIC;
    }
  }

  return Photon::GLOBAL;
}

Photon::Photon(R3Ray ray, RNRgb rgb, int photon_type, Photon *source=NULL)
  : ray(ray), source(source), rgb(rgb),
    photon_type(photon_type), path_type(GetPathType(photon_type, source))
{
}

void
Photon::Draw(double radius) const
{
  glColor3d(rgb.R(), rgb.G(), rgb.B());

  R3Point start = this->ray.Start();
  R3Vector dir = this->ray.Vector();

  R3Span(start, start + dir * radius).Draw();
}

void
Photon::DrawPath(double radius) const
{
  const Photon *photon = this;
  while (photon != NULL && photon->source != NULL) {
    R3Point start = photon->ray.Start();
    R3Point end = photon->source->ray.Start();
    RNRgb source_rgb = photon->source->rgb;

    glColor3d(source_rgb.R(), source_rgb.G(), source_rgb.B());
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

R3Vector
RotateToVector(R3Vector v, R3Vector original_norm, R3Vector new_norm)
{
  original_norm.Normalize();
  new_norm.Normalize();

  R3Vector axis = R3Vector(original_norm);
  axis.Cross(new_norm);
  axis.Normalize();

  RNScalar cos_angle = new_norm.Dot(original_norm);
  RNAngle angle = acos(cos_angle);

  v.Rotate(axis, angle);

  return v;
}

Photon *
PhotonFromDirLight(R3DirectionalLight *light, int scene_radius)
{
  R3Vector dir = light->Direction();
  dir.Normalize();

  RNScalar x, z;
  while (true) {
    x = 1.5 * scene_radius * (2.0 * Random() - 1.0);
    z = 1.5 * scene_radius * (2.0 * Random() - 1.0);
    if (x*x + z*z < 2 * scene_radius * scene_radius) break;
  }

  R3Vector pos = R3Vector(x,0,z);
  pos = RotateToVector(pos, R3Vector(0,1,0), light->Direction());
  pos -= scene_radius * 2 * dir;

  R3Ray ray = R3Ray(pos.Point(), dir);
  Photon *photon = new Photon(ray, light->Color(), Photon::LIGHT);

  return photon;
}
Photon *
PhotonFromPointLight(R3PointLight *light)
{
  R3Ray ray = R3Ray(light->Position(), RandomVectorUniform());
  Photon *photon = new Photon(ray, light->Color(), Photon::LIGHT);

  return photon;
}
Photon *
PhotonFromSpotLight(R3SpotLight *light)
{
  R3Ray ray = R3Ray(light->Position(), RandomVectorSpot(light));
  Photon *photon = new Photon(ray, light->Color(), Photon::LIGHT);

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

R3Ray
DiffuseBounce(R3Point pos, R3Vector norm) {
  R3Vector dir = RandomVectorInDir(norm);
  return R3Ray(pos, dir);
}

R3Ray
SpecularBounce(R3Vector source_dir, R3Point inters_pos, R3Vector inters_norm) {
  R3Vector dir = 2 * inters_norm * inters_norm.Dot(- source_dir) + source_dir;

  return R3Ray(inters_pos, dir);
}

R3Ray
TransmissionBounce(R3Vector source_dir,
    RNScalar index_of_refraction,
    R3Point inters_pos, R3Vector inters_norm) {

  R3Vector l = - source_dir;
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

  return R3Ray(inters_pos, dir);
}

RNScalar
MaxRGB(RNRgb c) {
  return fmax(fmax(c.R(), c.G()), c.B());
}

RNScalar
CoeffScaledByColor(RNRgb coeff, RNRgb pow) {
  return MaxRGB(pow * coeff) / MaxRGB(pow);
}

void
GetCoeffs(RNRgb rgb, const R3Brdf *brdf,
  RNRgb &diff, RNRgb &spec, RNRgb &trans,
  RNScalar &kd, RNScalar &ks, RNScalar &kt, RNScalar &k_total)
{
  diff = ScaleRGB(brdf->Diffuse());
  spec = ScaleRGB(brdf->Specular());
  trans = ScaleRGB(brdf->Transmission());

  // diff = brdf->Diffuse();
  // spec = brdf->Specular();
  // trans = brdf->Transmission();

  kd = CoeffScaledByColor(diff, rgb);
  ks = CoeffScaledByColor(spec, rgb);
  kt = CoeffScaledByColor(trans, rgb);

  k_total = fmax(1.0, kd + ks + kt);
}

RNBoolean
Bounce(R3Vector source_dir, RNRgb source_rgb,
       R3Point point, R3Vector normal,
       R3Material *material,
       R3Ray &new_ray, RNRgb &new_rgb, int &collision_type)
{
  const R3Brdf *brdf = material->Brdf();

  RNRgb diff, spec, trans;
  RNScalar kd, ks, kt, k_total;
  GetCoeffs(source_rgb, brdf, diff, spec, trans, kd, ks, kt, k_total);

  RNScalar r = Random() * k_total + 0.2;
  if (r < kd) {
    new_rgb = diff * source_rgb / kd; // scale inverse to probability
    collision_type = Photon::DIFFUSE;
    new_ray = DiffuseBounce(point, normal);
    return true;

  } else if (r < ks + kd) {
    new_rgb = spec * source_rgb / ks;
    collision_type = Photon::SPECULAR;
    new_ray = SpecularBounce(source_dir, point, normal);
    return true;

  } else if (r < kt + ks + kd) {
    RNScalar ir = brdf->IndexOfRefraction();
    new_rgb = trans * source_rgb / kt;
    collision_type = Photon::TRANSMISSION;
    new_ray = TransmissionBounce(source_dir, ir, point, normal);
    return true;

  }

  return false;
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
  R3Ray source_ray = source_photon->ray;
  R3Vector source_dir = source_ray.Vector();
  R3Ray ray = R3Ray(source_ray.Point(0.01), source_dir);

  RNBoolean intersected = scene->Intersects(ray,
    &node, &element, &shape, &point, &normal, &t);

  if (!intersected) {
    return NULL;
  }

  RNRgb new_rgb;
  int collision_type;
  R3Ray new_ray;
  RNBoolean bounced = Bounce(source_dir, source_photon->rgb,
      point, normal, element->Material(), new_ray, new_rgb, collision_type);

  if (bounced) {
    return new Photon(new_ray, new_rgb, collision_type, source_photon);
  } else {
    return NULL;
  }
}

////////////////////////////////////////////////////////////////////////
// Store the photons in a KDTree and render them!
////////////////////////////////////////////////////////////////////////

PhotonSample *
ProcessPhotonSample(Photon *photon) {
  if (photon->source == NULL) {
    return NULL;
  }

  R3Point position = photon->ray.Start();
  R3Point source = photon->source->ray.Start();
  R3Vector incident = position - source;
  incident.Normalize();

  PhotonSample *photon_sample = new PhotonSample;
  photon_sample->position = position;
  photon_sample->incident = incident;

  photon_sample->rgb = photon->rgb;

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

  int initial_photons = 1000000;
  int total_photons = initial_photons * 10;

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

RNRgb
SamplePhotonsAtPoint(R3Scene *scene, R3Point point, R3Vector normal)
{
  RNRgb pixel_color = RNRgb(0,0,0);

  // Load samples
  R3Kdtree<PhotonSample *> *global_samples = GetKdGlobalSamples(scene);
  R3Kdtree<PhotonSample *> *caustic_samples = GetKdCausticSamples(scene);

  RNArray<PhotonSample *> samples = RNArray<PhotonSample *>();

  double radius = scene->BBox().DiagonalRadius();
  global_samples->FindClosest(point, 0, radius * 0.05, 500, samples);
  caustic_samples->FindClosest(point, 0, radius * 0.01, 500, samples);

  RNScalar num_samples = samples.NEntries();
  for (int i = 0; i < num_samples; i ++) {
    PhotonSample *sample = samples.Kth(i);
    RNRgb sample_rgb = sample->rgb;

    // If we assume each of these flux samples is passing through the
    // intersection point, how much of that flux is reflected towards the viewer?

    R3Vector incident_dir = -sample->incident;
    R3Vector specular_dir = R3Vector(incident_dir);
    specular_dir.Rotate(normal, M_PI);

    RNScalar diffuse_scaling = fmax(0.0, incident_dir.Dot(normal));
    pixel_color += sample_rgb * diffuse_scaling / num_samples;
  }

  return ClampRGB(pixel_color);
}

RNRgb RenderRay(R3Scene *scene, R3Ray source_ray) {
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  // offset the ray a bit
  R3Vector source_dir = source_ray.Vector();
  R3Ray ray = R3Ray(source_ray.Point(0.01), source_dir);

  RNBoolean intersected = scene->Intersects(ray,
    &node, &element, &shape, &point, &normal, &t);

  if (!intersected) {
    return scene->Ambient();
  }

  R3Material *material = element->Material();
  const R3Brdf *brdf = material->Brdf();

  // get diffuse color at intersection
  RNRgb rgb = SamplePhotonsAtPoint(scene, point, normal);

  RNRgb diff, spec, trans;
  RNScalar kd, ks, kt, k_total;
  GetCoeffs(rgb, brdf, diff, spec, trans, kd, ks, kt, k_total);

  rgb = diff * rgb * kd / k_total;

  // do specular/transmission bounce
  // RNScalar r = Random() * fmax(1.0, ks + kt);
  // if (r < ks) {
  //   R3Ray new_ray = SpecularBounce(source_dir, point, normal);
  //   RNRgb new_rgb = RenderRay(scene, new_ray);
  //   rgb += spec * new_rgb * ks / k_total;
  //
  // } else if (r < ks + kt) {
  //   RNScalar ir = brdf->IndexOfRefraction();
  //   R3Ray new_ray = TransmissionBounce(source_dir, ir, point, normal);
  //   RNRgb new_rgb = RenderRay(scene, new_ray);
  //   rgb += trans * new_rgb * kt / k_total;
  // }

  return ClampRGB(rgb);
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
      RNRgb color = RenderRay(scene, ray);
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
