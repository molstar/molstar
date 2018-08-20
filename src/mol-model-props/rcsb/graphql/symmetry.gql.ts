 // workaround so the query gets found by the codegen
function gql (strs: TemplateStringsArray) { return strs.raw.join('') }

export default
gql`query AssemblySymmetry($pdbId: String!) {
    assemblies(pdbId: $pdbId) {
        assembly_id
        rcsb_assembly_symmetry {
            source
            symmetry_features {
                symmetry {
                    description
                    value
                }
                clusters {
                    members
                    avg_rmsd
                }
                stoichiometry {
                    description
                    value
                }
                symmetry_axes {
                    start
                    end
                    order
                }
                type
            }
        }
    }
}`