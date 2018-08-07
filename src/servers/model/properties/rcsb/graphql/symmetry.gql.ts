 // workaround so the query gets found by the codegen
function gql (strs: TemplateStringsArray) { return strs.raw.join('') }

export default
gql`query RcsbSymmetry($pdbId: String) {
    assemblies(pdbId: $pdbId) {
        assembly_id
        rcsb_annotation_symmetry {
            source
            symmetry_features {
                type
                clusters {
                    avg_rmsd
                    members
                }
                stoichiometry {
                    description
                    value
                }
                symmetry_axes {
                    order
                    start
                    end
                }
            }
        }
    }
}`