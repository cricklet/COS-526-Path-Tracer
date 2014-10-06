import os

for res in [200]:
 for light_photons in [100000]:
  for pixel_samples in [4, 16, 64]:
   for photon_samples in [64]:
    for sample_dist in [0.05]:
      input_name = 'cornell'
      output_name = '%s_%d_%d_%d_%f_%d' % (input_name, res, light_photons, photon_samples, sample_dist, pixel_samples)
      params = '-resolution %d %d -light_photons %d -photon_samples %d -sample_dist %f -pixel_samples %d'\
               % (res, res, light_photons, photon_samples, sample_dist, pixel_samples)
      cmd = './src/photonmap ./input/%s.scn -v %s ./output/%s.bmp' % (input_name, params, output_name)

      print cmd
      os.system(cmd)
      
