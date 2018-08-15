// /**
//  * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author David Sehnal <david.sehnal@gmail.com>
//  */

// export default [{
//     name: 'Residues connected to HEM',
//     value: `(sel.atom.is-connected-to
//   sel.atom.res
//   :target (sel.atom.res (= atom.label_comp_id HEM))
//   ;; default bond test allows only covalent bonds
//   :bond-test true
//   :disjunct true)`
// }, {
//     name: 'All C or N atoms in ALA residues',
//     value: `(sel.atom.atom-groups
//   :residue-test (= atom.auth_comp_id ALA)
//   :atom-test (set.has (set _C _N) atom.el))`
// }, {
//     name: 'Residues 130 to 180',
//     value: `(sel.atom.res (in-range atom.resno 130 180))`
// }, {
//     name: 'All residues within 5 ang from Fe atom',
//     value: `(sel.atom.include-surroundings
//   (sel.atom.atoms (= atom.el _Fe))
//   :radius 5
//   :as-whole-residues true)`
// }, {
//     name: 'Cluster LYS residues within 5 ang',
//     value: `(sel.atom.cluster
//   (sel.atom.res (= atom.label_comp_id LYS))
//   :max-distance 5)`
// }, {
//     name: 'Residues with max b-factor < 45',
//     value: `(sel.atom.pick sel.atom.res
//   :test (< (atom.set.max atom.B_iso_or_equiv) 45))`
// }, {
//   name: 'Atoms between 10 and 15 ang from Fe',
//   value: `(sel.atom.within sel.atom.atoms
//   :target (sel.atom.atoms (= atom.el _Fe))
//   :min-radius 10
//   :max-radius 15)`
// }, {
//   name: 'HEM and 2 layers of connected residues',
//   value: `(sel.atom.include-connected
//   (sel.atom.res (= atom.label_comp_id HEM))
//   ;; default bond test allows only covalent bonds
//   ;; another option is to use :bond-test true to allow any connection
//   :bond-test (bond.is metallic covalent)
//   :layer-count 2
//   :as-whole-residues true)`
// }, {
//   name: 'All rings',
//   value: `(sel.atom.rings)`
// }, {
//   name: 'CCCCN and CCNCN rings',
//   value: `(sel.atom.rings
//   (ringfp _C _N _C _N _C)
//   ;; the "rotation" of element symbols has no effect
//   ;; the following is the same as (ringfp _C _C _C _C _N)
//   (ringfp _C _C _C _N _C))`
// }, {
//   name: 'Sheets',
//   value: `(sel.atom.atom-groups
//   :residue-test (atom.sec-struct.is sheet)
//   :group-by (atom.key.sec-struct))`
// }, {
//   name: 'Helices formed by at least 30 residues',
//   value: `(sel.atom.pick
//   (sel.atom.atom-groups
//   :residue-test (atom.sec-struct.is helix)
//   :group-by atom.key.sec-struct)
//   :test (<= 30 (atom.set.count-query sel.atom.res)))`
// }, {
//   name: 'Modified residues',
//   value: `(sel.atom.res atom.is-modified)`
// }, {
//   name: 'Atoms participating in metallic coordination',
//   value: `(sel.atom.atoms
//   (> (atom.bond-count :flags (bond-flags metallic)) 0))`
// }, {
//   name: 'LYS and ALA residues that are between 2 and 5 ang apart',
//   value: `(sel.atom.dist-cluster
//   ;; upper triangular matrix are maximum distances of corresponding selections
//   ;; lower triangular matrix are minumum distances of corresponding selections
//   :matrix [
//     [0 5]
//     [2 0]
//   ]
//   :selections [
//     (sel.atom.res (= atom.resname LYS))
//     (sel.atom.res (= atom.resname ALA))
//   ])`
// }, {
//   name: 'Clusters of 3 LYS residues that are mutually no more than 5 ang apart',
//   value: `(sel.atom.dist-cluster
//   :matrix [[0 5 5] [0 0 5] [0 0 0]]
//   :selections [
//     (sel.atom.res (= atom.resname LYS))
//     (sel.atom.res (= atom.resname LYS))
//     (sel.atom.res (= atom.resname LYS))
//   ])`
// }]