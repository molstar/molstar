## Particles

- A target's atoms colored by element near the camera, blended in from its particle's color by their size on screen
  (the target representation's `structure.siteColor`: on by default, carbons keeping the particle's color, a radius of
  4 to 12 px; WebGL and WebGPU, the sphere raster too). For targets of known elements: a structure target, an mmCIF
  assembly's sites, or a particle pack's of version 5. A pack's by site type, no data a site: the yeast cell's 7 types;
  a structure target's a color an atom.
- Particle packs version 5: a target's sites keep their elements (`siteElements`, a type a distinct radius and
  element; `ParticleSites.elements`, `getElement`, `hasElements`). Versions 2 to 4 still read, without elements: packed
  again to have them.
