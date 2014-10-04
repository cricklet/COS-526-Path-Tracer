// Source file for photonmap renderer
// This implementation is a simple raycaster.
// Replace it with your own code.



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Graphics/R3Graphics.h"

////////////////////////////////////////////////////////////////////////
// Helpers.
////////////////////////////////////////////////////////////////////////

RNScalar Random() {
  return static_cast <RNScalar> (rand())
       / static_cast <RNScalar> (RAND_MAX);
}

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

RNScalar
MaxRGB(RNRgb c) {
  return fmax(fmax(c.R(), c.G()), c.B());
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

void
GetBRDFProbabilities(RNRgb source_rgb, const R3Brdf *brdf,
  RNRgb& diff, RNRgb& spec, RNRgb& trans,
  RNScalar& kd, RNScalar& ks, RNScalar& kt)
{
  diff  = brdf->Diffuse();
  spec  = brdf->Specular();
  trans = brdf->Transmission();
  RNScalar total = MaxRGB(diff) + MaxRGB(spec) + MaxRGB(trans);

  if (total > 0.8) {
    diff  /= total;
    spec  /= total;
    trans /= total;
    diff  *= 0.8;
    spec  *= 0.8;
    trans *= 0.8;
  }

  if (MaxRGB(source_rgb) == 0) {
    kd = 0;
    ks = 0;
    kt = 0;
  } else {
    kd = MaxRGB(source_rgb * diff);//  / MaxRGB(source_rgb);
    ks = MaxRGB(source_rgb * spec);//  / MaxRGB(source_rgb);
    kt = MaxRGB(source_rgb * trans);// / MaxRGB(source_rgb);
  }
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

class PhotonPath {
public:
  enum { LIGHT, DIFFUSE, SPECULAR, TRANSMISSION };

  PhotonPath(R3Point start, R3Point end, RNRgb rgb,
             int start_type, PhotonPath *next);

  void Draw(double radius, RNBoolean draw_entire_path) const;

  R3Vector Incident() const;

  const R3Point start;
  const R3Point end;
  const RNRgb rgb;
  const int start_type;
  const PhotonPath *next;
};

PhotonPath::PhotonPath(R3Point start, R3Point end,
  RNRgb rgb, int start_type,
  PhotonPath *next)
  : start(start), end(end),
    rgb(rgb), start_type(start_type),
    next(next)
{
}

R3Vector
PhotonPath::Incident() const
{
  R3Vector incident = end - start;
  incident.Normalize();
  return incident;
}

void
PhotonPath::Draw(double radius, RNBoolean draw_entire_path) const
{
  glColor3d(rgb.R(), rgb.G(), rgb.B());

  R3Vector incident = this->Incident();
  R3Span(end, end - incident * radius).Draw();

  if (draw_entire_path) {
    R3Span(start, end).Draw();
  }
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

RNBoolean
Bounce(R3Vector source_dir, RNRgb source_rgb,
       R3Point point, R3Vector normal,
       R3Material *material,
       R3Ray &new_ray, RNRgb &new_rgb, int &collision_type)
{
  const R3Brdf *brdf = material->Brdf();

  RNRgb diff, spec, trans;
  RNScalar kd, ks, kt;
  GetBRDFProbabilities(source_rgb, brdf, diff, spec, trans, kd, ks, kt);

  RNScalar r = Random();
  if (r < kd) {
    new_rgb = source_rgb * diff / kd; // scale inverse to probability
    collision_type = PhotonPath::DIFFUSE;
    new_ray = DiffuseBounce(point, normal);
    return true;

  } else if (r < ks + kd) {
    new_rgb = source_rgb * spec / ks;
    collision_type = PhotonPath::SPECULAR;
    new_ray = SpecularBounce(source_dir, point, normal);
    return true;

  } else if (r < kt + ks + kd) {
    RNScalar ir = brdf->IndexOfRefraction();
    new_rgb = source_rgb * trans / kt;
    collision_type = PhotonPath::TRANSMISSION;
    new_ray = TransmissionBounce(source_dir, ir, point, normal);
    return true;

  }

  return false;
}

////////////////////////////////////////////////////////////////////////
// Randomly sample light sources based on their intensity.
// Generate photons from light sources.
////////////////////////////////////////////////////////////////////////

R3Ray
RayFromDirLight(R3DirectionalLight *light, int scene_radius)
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

PhotonPath *
CreatePhotonPath(R3Scene *scene, R3Ray start_ray, RNRgb start_color, int start_type)
{
  // intersection variables
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point intersection;
  R3Vector normal;
  RNScalar t;

  // offset the ray a bit
  R3Ray ray = R3Ray(start_ray.Point(0.01), start_ray.Vector());

  RNBoolean intersected = scene->Intersects(ray,
    &node, &element, &shape, &intersection, &normal, &t);

  // don't create a photon path if there's no intersection
  if (!intersected) {
    return NULL;
  }

  // bounce at the intersection
  R3Ray new_ray;
  RNRgb new_rgb;
  int new_type;
  RNBoolean bounced = Bounce(
    start_ray.Vector(), start_color, intersection, normal,
    element->Material(), new_ray, new_rgb, new_type);

  // create a new photon path at the intersection
  PhotonPath *new_path = NULL;
  if (bounced) {
    new_path = CreatePhotonPath(scene, new_ray, new_rgb, new_type);
  }

  // create path from start to intersection, with start's color/type
  PhotonPath *path = new PhotonPath(start_ray.Start(), intersection,
    start_color, start_type, new_path);

  return path;
}

RNArray<PhotonPath *> *
CreatePhotonPathsFromLights(R3Scene *scene, int num)
{
  printf("Generating photon paths from lights.\n");
  RNArray<PhotonPath *> *paths = new RNArray<PhotonPath *>;

  RNScalar total_intensity = TotalLightIntensity(scene);
  printf("Total light intensity %f.\n", total_intensity);

  double radius = scene->BBox().DiagonalRadius();
  for (int i = 0; i < num; i ++) {
    R3Light *light = RandomLight(scene, total_intensity);

    int light_class = light->ClassID();

    R3Ray ray;

    if (light_class == R3DirectionalLight::CLASS_ID()) {
      ray = RayFromDirLight((R3DirectionalLight *) light, radius);
    }
    else if (light_class == R3PointLight::CLASS_ID()) {
      ray = RayFromPointLight((R3PointLight *) light);
    }
    else if (light_class == R3SpotLight::CLASS_ID()) {
      ray = RayFromSpotLight((R3SpotLight *) light);
    }

    PhotonPath *path = CreatePhotonPath(scene, ray,
      light->Color(), PhotonPath::LIGHT);

    paths->Insert(path);
  }

  return paths;
}

////////////////////////////////////////////////////////////////////////
// Store the photons in a KDTree and render them!
////////////////////////////////////////////////////////////////////////

PhotonSample *
CreatePhotonSample(const PhotonPath *path) {
  R3Vector incident = path->Incident();

  PhotonSample *photon_sample = new PhotonSample;
  photon_sample->position = path->end;
  photon_sample->incident = incident;
  photon_sample->rgb = path->rgb;

  return photon_sample;
}

////////////////////////////////////////////////////////////////////////
// Function for generating + caching photons
// (both path tracing + rendering)
////////////////////////////////////////////////////////////////////////

static RNArray<PhotonPath *> *cached_paths = NULL;

static RNArray<PhotonSample *> *cached_global_samples = NULL;
static R3Kdtree<PhotonSample *> *cached_kd_global_samples = NULL;

static RNArray<PhotonSample *> *cached_caustic_samples = NULL;
static R3Kdtree<PhotonSample *> *cached_kd_caustic_samples = NULL;

void
CachePhotonPaths(R3Scene *scene) {
  if (cached_paths  != NULL) {
    return;
  }

  RNArray<PhotonPath *> *paths = CreatePhotonPathsFromLights(scene, 100000);
  printf("Created initial %d photons\n", paths->NEntries());

  cached_paths = paths;
}

RNArray<PhotonPath *> *
GetPhotonPaths(R3Scene *scene) {
  CachePhotonPaths(scene);
  return cached_paths;
}

void
CacheSamples(R3Scene *scene) {
  if (cached_global_samples     != NULL
      && cached_caustic_samples != NULL) {
    return;
  }
  printf("Generating photon samples\n");

  RNArray<PhotonSample *> *global_samples = new RNArray<PhotonSample *>;
  RNArray<PhotonSample *> *caustic_samples = new RNArray<PhotonSample *>;

  RNArray<PhotonPath *> *paths = GetPhotonPaths(scene);
  for (int i = 0; i < paths->NEntries(); i ++) {
    const PhotonPath *path = paths->Kth(i);
    while (path != NULL) {
      if (path->start_type == PhotonPath::TRANSMISSION
          || path->start_type == PhotonPath::SPECULAR) {
        caustic_samples->Insert(CreatePhotonSample(path));
      } else {
        global_samples->Insert(CreatePhotonSample(path));
      }
      path = path->next;
    }
  }

  R3Kdtree<PhotonSample *> *kd_global_samples =
    new R3Kdtree<PhotonSample *>(*global_samples, GetPhotonSamplePosition);
  R3Kdtree<PhotonSample *> *kd_caustic_samples =
    new R3Kdtree<PhotonSample *>(*caustic_samples, GetPhotonSamplePosition);

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
DrawPhotonPaths(R3Scene *scene)
{
  double radius = scene->BBox().DiagonalRadius();

  RNArray<PhotonPath *> *paths = GetPhotonPaths(scene);
  for (int i = 0; i < paths->NEntries(); i ++) {
    const PhotonPath *path = paths->Kth(i);
    RNBoolean draw_entire_path = Random() > 0.999;

    while (path != NULL) {
      path->Draw(radius * 0.01, draw_entire_path);
      path = path->next;
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Function to render image with photon mapping
////////////////////////////////////////////////////////////////////////

RNRgb
SamplePhotonsAtPoint(R3Scene *scene,
  R3Point point, R3Vector normal, const R3Brdf *brdf)
{
  RNRgb pixel_color = RNRgb(0,0,0);

  // Load samples
  R3Kdtree<PhotonSample *> *global_samples = GetKdGlobalSamples(scene);
  R3Kdtree<PhotonSample *> *caustic_samples = GetKdCausticSamples(scene);

  RNArray<PhotonSample *> samples = RNArray<PhotonSample *>();

  double radius = scene->BBox().DiagonalRadius();
  global_samples->NNodes();
  global_samples->FindClosest(point, 0, radius * 0.1, 100, samples);
  caustic_samples->FindClosest(point, 0, radius * 0.05, 100, samples);

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

    RNRgb diffuse_coeff = ScaleRGB(brdf->Diffuse());

    pixel_color += sample_rgb * diffuse_scaling * diffuse_coeff / 100.0;
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
  RNRgb rgb = SamplePhotonsAtPoint(scene, point, normal, brdf);

  // do specular & transmission bounce
  RNRgb diff, spec, trans;
  RNScalar kd, ks, kt;
  GetBRDFProbabilities(RNRgb(1,1,1), brdf, diff, spec, trans, kd, ks, kt);

  if (kt > 0) {
    printf("asdfasdf\n");
  }

  RNScalar r = Random();

  if (r < ks) {
    R3Ray spec_ray = SpecularBounce(source_dir, point, normal);
    RNRgb spec_rgb = RenderRay(scene, spec_ray);
    rgb += spec * spec_rgb;
  } else if (r < ks + kt) {
    RNScalar ir = brdf->IndexOfRefraction();
    R3Ray trans_ray = TransmissionBounce(source_dir, ir, point, normal);
    RNRgb trans_rgb = RenderRay(scene, trans_ray);
    rgb += spec * trans_rgb;
  }

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
      RNRgb color = RNRgb(0,0,0);

      int num_samples = 8;
      for (int k = 0; k < num_samples; k ++) {
        color += RenderRay(scene, ray) / num_samples;
      }
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
