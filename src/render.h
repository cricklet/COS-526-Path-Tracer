// Include file for the photon map render code

R2Image *RenderImage(R3Scene *scene, int width, int height, int print_verbose);
void DrawPhotonPaths(R3Scene *scene);
void DrawCausticSamples(R3Scene *scene);
void DrawGlobalSamples(R3Scene *scene);

void SetNumLightPhotons(int n);
void SetNumPixelSamples(int n);
void SetNumPhotonSamples(int n);
void SetPhotonSampleDist(double d);
