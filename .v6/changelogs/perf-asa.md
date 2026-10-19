## Model: faster accessible surface area

- Shrake-Rupley accessible surface area no longer allocates per atom (typed neighbour buffers, no sort, last occluder
  first); the areas are unchanged. 186k atoms: 2.9 s → 1.7 s, 5.6 GB → 11 MB allocated; 7qpd (11.5k atoms): 2× faster
  (#1986)
