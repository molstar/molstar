## Graphics: many small units

- [Breaking] `Sphere3D.extremaRadius`: a set of at most 98 points with radii keeps its points as bounding-sphere
  extrema with their largest radius, a sphere keeps no duplicate extrema, and `Sphere3D.expand` grows that radius
  instead of moving the extrema; readers of `extrema` add it (`Sphere3D.getExtremaRadius`)
- 100k small units as merged spacefill with the per-unit visual: built in 6.7 s with 2.1 GB of heap (did not fit
  before); 10k units: 1.54 → 0.45 s, 43.8 → 23.0 KB a unit; a 94k-unit NPC capsid: 17.5 → 9.1 s, 4.6 → 2.8 GB
- Many small units: layers without a grid share their placeholder texture and grid values, a single instance gets a
  one-cell grid without the grid lookups, merge keys are shared by equal members; 100k units 6.7 → 3.8 s
