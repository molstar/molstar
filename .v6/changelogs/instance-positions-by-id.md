## Graphics

- Fix per-vertex colors and sizes of instanced geometry (e.g. a volume color theme on an assembly's surface), and the
  instance colors and sizes of geometry exports (glTF, OBJ, STL, USDZ), being placed by another instance's transform
  where the instance grid reorders the instances: the position iterators and the exports read an instance's transform
  at its row (`getInstanceRows`, `getInstanceMatricesById`)
