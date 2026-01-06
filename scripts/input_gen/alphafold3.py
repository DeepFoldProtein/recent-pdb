#!/usr/bin/env python3
"""
AlphaFold3 Input Generator
==========================
Generates JSON input format for AlphaFold3 inference.

AlphaFold3 expects a specific JSON schema with:
- sequences: List of chain sequences with entity information
- msas: MSA data per entity
- templates: Template structures (optional)
"""

import argparse

import json
from pathlib import Path

import yaml

from .base import (
    ModelInputGenerator,
    TargetFeatures,
    load_target_features,
    logger,
)


class AlphaFold3InputGenerator(ModelInputGenerator):
    """Generate AlphaFold3 JSON input format."""

    name = "AlphaFold3"
    input_format = "json"

    def generate(
        self,
        target_features: TargetFeatures,
        output_dir: Path,
    ) -> Path:
        """Generate AlphaFold3 input JSON."""
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{target_features.target_id}.json"

        # Build input structure
        af3_input = {
            "name": target_features.target_id,
            "dialect": "alphafold3",
            "version": 1,
            "modelSeeds": [42],  # Default seed
            "sequences": [],
        }

        # Group chains by sequence hash (for homo-oligomers)
        hash_to_chains: dict[str, list] = {}
        for chain in target_features.chains:
            if chain.seq_hash not in hash_to_chains:
                hash_to_chains[chain.seq_hash] = []
            hash_to_chains[chain.seq_hash].append(chain)

        # Build sequence entries
        for seq_hash, chains in hash_to_chains.items():
            first_chain = chains[0]

            # Load and combine MSAs
            msa_data = []
            seen_seqs = set()
            templates = []

            # Add query sequence first
            if first_chain.sequence:
                seen_seqs.add(first_chain.sequence)
                msa_data.append({"sequence": first_chain.sequence})

            if first_chain.msa_paths:
                for path in first_chain.msa_paths:
                    # Check if this is a template database (pdb100)
                    if "pdb100" in str(path) and not templates:
                        templates = self._process_templates(path, first_chain.sequence)
                        continue

                    # Otherwise treat as MSA
                    if "pdb100" not in str(path):
                        msa_seqs = self.load_msa(path)
                        for _, seq in msa_seqs:
                            # Standardize A3M -> Aligned FASTA (remove insertions)
                            # Keep uppercase (matches) and dashes (deletions)
                            seq_clean = "".join(
                                c for c in seq if c.isupper() or c == "-"
                            )
                            if seq_clean not in seen_seqs:
                                seen_seqs.add(seq_clean)
                                msa_data.append({"sequence": seq_clean})

            # Limit MSA depth
            msa_data = msa_data[:2048]

            # Add MSA if available
            if msa_data:
                unpaired_msa = "\n".join(
                    f">{i}\n{m['sequence']}" for i, m in enumerate(msa_data)
                )
            else:
                unpaired_msa = f">{first_chain.sequence}\n{first_chain.sequence}"

            sequence_entry = {
                "protein": {
                    "id": [c.chain_id for c in chains],
                    "sequence": first_chain.sequence,
                    "unpairedMsa": unpaired_msa,
                    "pairedMsa": "",
                    "templates": templates,
                },
            }

            af3_input["sequences"].append(sequence_entry)

        # Write output
        with open(output_path, "w") as f:
            json.dump(af3_input, f, indent=2)

        logger.info(f"Generated AlphaFold3 input: {output_path}")
        return output_path

    def _process_templates(self, a3m_path: Path, query_sequence: str) -> list[dict]:
        """Parse PDB100 A3M and generate AF3 template objects."""
        # query_sequence unused but kept for API potential
        _ = query_sequence
        templates = []
        try:
            import gemmi
        except ImportError:
            logger.warning("gemmi not installed, skipping templates")
            return []

        import os

        logger.info(f"Processing templates from: {a3m_path}")
        pdb_master = Path(os.path.expandvars(self.config["paths"]["pdb_master"]))

        # Extract matches from A3M
        # Limit to top 4 templates to avoid huge input
        matches = []
        with open(a3m_path) as f:
            lines = f.readlines()

        # A3M format: Query first, then hits.
        # >Query
        # SEQ
        # >Hit
        # SEQ
        if len(lines) < 4:
            return []

        # Parse query from A3M check?
        # Assuming standard mmseqs output

        # Process hits
        seen_pdbs = set()

        # Skip query (first 2 lines if 2-line fasta-ish)
        # Actually load_msa logic helps?
        # Let's iterate manually to access headers

        # Skip first entry (Query)
        # Identify by being first.
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith(">"):
                header = line[1:]
                seq_lines = []
                i += 1
                while i < len(lines) and not lines[i].startswith(">"):
                    seq_lines.append(lines[i].strip())
                    i += 1
                seq = "".join(seq_lines)

                seq = "".join(seq_lines)
                matches.append((header, seq))
            else:
                i += 1

        # Drop first match (Query)
        if matches:
            matches.pop(0)

        logger.info(f"Found {len(matches)} potential templates in A3M")

        for header, a3m_seq in matches:
            if len(templates) >= 4:
                break

            # Parse Header: >pdb_id_Chain ...
            # e.g. >2h35_C ...
            parts = header.split()[0].split("_")
            if len(parts) < 2:
                continue

            pdb_id = parts[0].lower()
            chain_id = parts[1]

            if pdb_id in seen_pdbs:
                continue

            # Find CIF
            cif_path = pdb_master / f"{pdb_id[1:3]}/{pdb_id}.cif.gz"
            if not cif_path.exists():
                cif_path = pdb_master / f"{pdb_id[1:3]}/{pdb_id}.cif"
            if not cif_path.exists():
                # Try typical divisions?
                # Assuming simple structure or flat?
                # config usually points to divided/mmCIF/
                # Try finding via glob if needed, but slow.
                # Let's assume standard divided layout (middle hash)
                cif_path = pdb_master / f"{pdb_id[1:3]}/{pdb_id}.cif.gz"

            if not cif_path.exists():
                logger.warning(f"Template CIF not found: {cif_path}")
                continue

            try:
                # Read CIF
                block = gemmi.cif.read(str(cif_path)).sole_block()

                # Verify chain exists
                # Extract sequence from CIF for mapping
                # We need to map A3M columns to indices

                # Simple mapping:
                # 1. Un-align A3M sequence (hit) -> pure sequence
                # 2. Find this sequence in the Structure

                hit_seq_ungapped = a3m_seq.replace("-", "").upper()
                # Remove insertions (lowercase)? mmseqs A3M lowercase = insertion relative to query.
                # But it IS present in the template. So keep it.
                hit_seq_ungapped = a3m_seq.replace("-", "").upper()

                # Get Entities
                # Use our sync_pdb logic or gemmi helper
                # We need the sequence of 'chain_id'

                # Find entity_id for chain_id
                # _struct_asym.id == chain_id -> entity_id
                target_entity_id = None
                asym = block.get_mmcif_category("_struct_asym.")
                if "id" in asym and "entity_id" in asym:
                    for aid, eid in zip(asym["id"], asym["entity_id"]):
                        if aid == chain_id:
                            target_entity_id = eid
                            break

                if not target_entity_id:
                    continue

                # Get sequence
                poly_seq = block.get_mmcif_category("_entity_poly_seq.")
                # ... extract seq ...
                # To save code space, let's trust gemmi's polymer handling if possible?
                # Actually, parsing _entity_poly_seq is robust.

                real_seq_list = []
                if "entity_id" in poly_seq and "mon_id" in poly_seq:
                    for eid, mon in zip(poly_seq["entity_id"], poly_seq["mon_id"]):
                        if eid == target_entity_id:
                            code = gemmi.one_letter_code([mon])
                            real_seq_list.append(code if code else "X")
                real_seq = "".join(real_seq_list)

                # Find offset
                start_idx = real_seq.find(hit_seq_ungapped)
                if start_idx == -1:
                    # alignment might be partial or have mutations
                    # Fallback or skip? Skip for high precision
                    continue

                # Build indices
                query_indices = []
                template_indices = []

                q_idx = 0  # Index in query_sequence
                t_idx = start_idx  # Index in real_seq

                # A3M string iteration
                # Query sequence (from arg) corresponds to A3M columns (uppercase/dash)
                # Lowercase are insertions in T (not in Q).

                # Wait, query_sequence passed in is the FULL query.
                # A3M might be local.
                # Is the First Sequence in A3M the FULL query?
                # Usually mmseqs result2msa outputs valid A3M where first seq is full query?
                # Let's assume yes.

                for char in a3m_seq:
                    if char == "-":
                        # Deletion in T relative to Q
                        # Q advances, T stays (gap in T)
                        q_idx += 1
                    elif char.isupper():
                        # Match/Mismatch
                        # Exists in both (aligned)
                        query_indices.append(q_idx)
                        template_indices.append(t_idx)
                        q_idx += 1
                        t_idx += 1
                    elif char.islower():
                        # Insertion in T sequence
                        # Not in Q
                        t_idx += 1
                        # q_idx stays
                    else:
                        # Unknown? skip
                        pass

                # Filter structure to single chain using gemmi
                # AF3 requires template CIF to contain ONLY the relevant chain
                try:
                    # 1. Read original to get date
                    orig_doc = gemmi.cif.read(str(cif_path))
                    orig_block = orig_doc.sole_block()
                    revision_date = orig_block.find_value(
                        "_pdbx_audit_revision_history.revision_date"
                    )
                    if not revision_date:
                        revision_date = "1970-01-01"

                    # 2. Read structure and select target chain
                    st = gemmi.read_structure(str(cif_path))

                    # Create a new structure for the single chain
                    new_st = gemmi.Structure()
                    new_model = gemmi.Model("1")
                    new_st.add_model(new_model, pos=0)

                    # Find the target chain in the first model of the original structure
                    found_chain = None
                    if len(st) > 0:
                        source_model = st[0]
                        for chain in source_model:
                            # Check name (label)
                            if chain.name == chain_id:
                                found_chain = chain
                                break

                    if found_chain:
                        try:
                            # Verify atoms before adding (debug)
                            logger.debug(
                                f"Adding chain {found_chain.name} with {found_chain.count_atom_sites()} atoms"
                            )
                            new_model.add_chain(found_chain)
                        except Exception as exc:
                            logger.warning(f"Direct chain add failed, cloning: {exc}")
                            new_model.add_chain(found_chain.clone())

                    # Add model to structure AFTER adding chain (crucial for atom preservation)
                    new_st.add_model(new_model, pos=0)

                    # Ensure entities are defined for AF3/gemmi to recognize polymer chains
                    new_st.setup_entities()

                    # Fix for sequence alignment (IndexError in AF3):
                    # 1. Inject full sequence from A3M.
                    # 2. Align atoms to this full sequence to fix label_seq_id.
                    try:
                        # hit_seq_ungapped comes from A3M
                        full_seq_str = hit_seq_ungapped

                        # Get structure's current sequence (from atoms)
                        res_list = list(new_model[0])
                        struc_seq_str = "".join(
                            [
                                gemmi.find_tabulated_residue(r.name).one_letter_code
                                for r in res_list
                            ]
                        )

                        # Python based alignment (struc -> full)
                        msg_debug = f"Aligning Struct ({len(struc_seq_str)}) to Full ({len(full_seq_str)})"
                        logger.debug(msg_debug)

                        struc_idx = 0
                        full_idx = 0
                        mapping = {}  # struc_res_idx -> full_seq_idx (0-based)

                        while struc_idx < len(struc_seq_str) and full_idx < len(
                            full_seq_str
                        ):
                            if struc_seq_str[struc_idx] == full_seq_str[full_idx]:
                                mapping[struc_idx] = full_idx
                                struc_idx += 1
                                full_idx += 1
                            else:
                                # Lookahead for match
                                found = False
                                for offset in range(1, 100):  # Lookahead limit
                                    if (
                                        full_idx + offset < len(full_seq_str)
                                        and struc_seq_str[struc_idx]
                                        == full_seq_str[full_idx + offset]
                                    ):
                                        full_idx += offset
                                        mapping[struc_idx] = full_idx
                                        struc_idx += 1
                                        full_idx += 1
                                        found = True
                                        break
                                if not found:
                                    # Forced mismatch consuming both
                                    mapping[struc_idx] = full_idx
                                    struc_idx += 1
                                    full_idx += 1

                        # Apply renumbering
                        for r_idx, seq_idx in mapping.items():
                            res = res_list[r_idx]
                            res.label_seq = seq_idx + 1  # 1-based index

                        # Now inject full sequence into entity
                        full_seq_codes = gemmi.expand_one_letter_sequence(
                            full_seq_str, gemmi.ResidueKind.AA
                        )
                        if len(new_st.entities) > 0:
                            new_st.entities[0].full_sequence = full_seq_codes

                        logger.debug("Renumbering and sequence injection complete.")

                    except Exception as e:
                        logger.warning(
                            f"Failed to isolate full sequence into entity: {e}"
                        )

                    if "found_chain" not in locals() or not found_chain:
                        logger.warning(f"Could not find chain {chain_id} in {cif_path}")
                        continue

                    # 3. Create new CIF and inject date
                    new_doc = new_st.make_mmcif_document()
                    new_block = new_doc.sole_block()

                    # Add revision history loop
                    loop = new_block.init_loop(
                        "_pdbx_audit_revision_history.", ["revision_date"]
                    )
                    loop.add_row([revision_date])

                    # MANUALLY GENERATE _pdbx_poly_seq_scheme
                    # Required by AF3 to link atoms to the entity sequence.
                    # Columns: asym_id, entity_id, seq_id, mon_id, pdb_seq_num, auth_seq_num, pdb_mon_id, auth_mon_id, pdb_strand_id, pdb_ins_code
                    try:
                        scheme_loop = new_block.init_loop(
                            "_pdbx_poly_seq_scheme.",
                            [
                                "asym_id",
                                "entity_id",
                                "seq_id",
                                "mon_id",
                                "pdb_seq_num",
                                "auth_seq_num",
                                "pdb_mon_id",
                                "auth_mon_id",
                                "pdb_strand_id",
                                "pdb_ins_code",
                            ],
                        )

                        # Create inverse mapping: full_seq_idx (0-based) -> res (gemmi.Residue)
                        # We need to re-derive this or capture it from the alignment block above.
                        # Since we renumbered label_seq to be 1-based seq_idx + 1, we can use that?
                        # Yes, res.label_seq corresponds to full sequence index + 1.

                        full_to_res = {}
                        res_list = list(
                            new_model[0]
                        )  # Get chains from model? No, new_model has residues?
                        # new_model structure: Model -> Chain -> Residue
                        # We added one chain.
                        chain_obj = new_model[0]
                        for r in chain_obj:
                            # r.label_seq should be set by our alignment logic
                            if r.label_seq:
                                full_to_res[r.label_seq - 1] = r

                        full_seq_codes = new_st.entities[0].full_sequence
                        ent_id = new_st.entities[0].name  # Usually "1"
                        asym_id = chain_obj.name  # Label ID? Or Name?
                        # In gemmi Chain.name is label_asym_id usually if loaded from CIF?
                        # But we created it or cloned it.
                        # gemmi.Structure.setup_entities assigns names.
                        # Let's use chain_obj.name for asym_id.

                        for i, mon in enumerate(full_seq_codes):
                            row = []
                            row.append(asym_id)
                            row.append(ent_id)
                            row.append(str(i + 1))  # seq_id
                            row.append(mon)  # mon_id

                            if i in full_to_res:
                                r = full_to_res[i]
                                row.append(
                                    str(r.seqid.num)
                                )  # pdb_seq_num - use original PDB numbering
                                row.append(
                                    str(r.seqid.num)
                                )  # auth_seq_num (same as pdb_seq_num)
                                row.append(r.name)  # pdb_mon_id
                                row.append(r.name)  # auth_mon_id
                                row.append(asym_id)  # pdb_strand_id
                                row.append(
                                    "?"
                                    if not r.seqid.icode or r.seqid.icode == " "
                                    else r.seqid.icode
                                )  # pdb_ins_code - use ? to match atoms
                            else:
                                # Start of gap
                                row.append("?")
                                row.append("?")
                                row.append("?")
                                row.append("?")
                                row.append("?")
                                row.append("?")

                            scheme_loop.add_row(row)

                        logger.debug(
                            f"Generated _pdbx_poly_seq_scheme with {len(full_seq_codes)} rows"
                        )

                    except Exception as e:
                        logger.warning(
                            f"[{pdb_id}] Failed to generate _pdbx_poly_seq_scheme: {e}"
                        )

                    cif_content = new_doc.as_string()

                    # Inject _atom_site.pdbx_PDB_model_num (required by AF3)
                    # Robust string patching moved AFTER as_string to prevent overwrite
                    import re

                    if "_atom_site.pdbx_PDB_model_num" not in cif_content:
                        # Insert tag after group_PDB
                        cif_content = cif_content.replace(
                            "_atom_site.group_PDB",
                            "_atom_site.group_PDB\n_atom_site.pdbx_PDB_model_num",
                        )
                        # Insert value '1' after ATOM/HETATM via regex
                        # Matches start of line ATOM/HETATM followed by space (auth seq id etc)
                        cif_content = re.sub(
                            r"^(ATOM|HETATM) ",
                            r"\1 1 ",
                            cif_content,
                            flags=re.MULTILINE,
                        )

                except Exception as e:
                    logger.warning(f"Error filtering CIF {cif_path}: {e}")
                    continue

                templates.append(
                    {
                        "mmcif": cif_content,
                        "queryIndices": query_indices,
                        "templateIndices": template_indices,
                    }
                )
                seen_pdbs.add(pdb_id)

            except Exception as e:
                logger.warning(f"Error processing template {pdb_id}: {e}")
                continue

        return templates


def main():
    parser = argparse.ArgumentParser(description="Generate AlphaFold3 input JSON")
    parser.add_argument("--target-id", required=True, help="Target PDB ID")
    parser.add_argument(
        "--cache-dir", type=Path, required=True, help="MSA cache directory"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Output directory"
    )
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))

    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    target_lists_dir = Path(config["paths"]["target_lists"])
    target_features = load_target_features(
        args.target_id, args.cache_dir, target_lists_dir, config
    )

    if target_features is None:
        return 1

    if target_features is None:
        return 1

    generator = AlphaFold3InputGenerator(config)
    output_path = generator.generate(target_features, args.output_dir / args.target_id)

    if generator.validate_output(output_path):
        logger.info("Input generation successful")
        return 0
    else:
        logger.error("Input generation failed")
        return 1


if __name__ == "__main__":
    exit(main())
