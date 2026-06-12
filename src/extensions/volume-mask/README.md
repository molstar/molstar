# Volume Mask Creator

Interactive 3D voxel mask generator for Mol*. Draw polygons from any camera angle to select regions of a cryo-EM volume, optionally combined with atomic structure proximity. Exports soft-edged masks as MRC/CCP4.

## Quick start (example app)

```bash
npm run dev -- -e volume-mask        # watch + build
http-server -p 1338 -g               # serve
# open http://localhost:1338/build/examples/volume-mask/
```
