# SteinhardtOrderParameters

`SteinhardtOrderParameters` computes per-atom bond-orientational order
parameters q_l, Wigner-3j invariants w_l and the normalized invariants
w_l_hat for a LAMMPS dump frame.

The implementation is a faithful port of LAMMPS `compute orientorder/atom`
(Thompson & Kohlmeyer), which in turn implements the order parameters
introduced by Steinhardt, Nelson & Ronchetti (1983) with the w_l_hat
normalization popularised by Lechner & Dellago (2008). q_6 and w_6_hat are
the canonical liquid/solid and crystallinity/nucleation discriminators.

## References

- P. J. Steinhardt, D. R. Nelson, M. Ronchetti, "Bond-orientational order
  in liquids and glasses", Phys. Rev. B 28, 784 (1983).
- W. Lechner, C. Dellago, "Accurate determination of crystal structures
  based on averaged local bond order parameters", J. Chem. Phys. 129,
  114707 (2008).
- LAMMPS `compute orientorder/atom` documentation
  (https://docs.lammps.org/compute_orientorder_atom.html).

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- SteinhardtOrderParameters
```

## CLI

Usage:

```bash
steinhardt-order-parameters <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--qlist <csv>` | No | Comma-separated list of orders l to compute. | `4,6,8,10,12` |
| `--nnn <int>` | No | Fixed number of nearest neighbors. `0` means use `--cutoff`. | `12` |
| `--cutoff <float>` | No | Cutoff radius (only used when `--nnn 0`). | `0.0` |
| `--wl <bool>` | No | Also compute w_l (Wigner-3j invariant). | `false` |
| `--wlHat <bool>` | No | Also compute w_l_hat = w_l / q_l^3 (dimensionless). | `false` |
| `--components <l>` | No | Emit the 2l+1 complex Q_l^m/|Q_l| components for a single l (disabled if < 0). | `-1` |
| `--onlySelected` | No | Only classify atoms flagged by a `Selection` column. | `false` |
| `--help` | No | Print CLI help. | |

> `--threads` is handled by the framework (OneTBB global_control), not by the plugin itself.

## Outputs

Two Parquet files are written next to `<output_base>`:

- `<output_base>_steinhardt.parquet` — summary (params + liquid/interface/crystal
  counts + a `structures` listing of structure_id/structure_name/atom_count).
- `<output_base>_atoms.parquet` — canonical D-009 per-atom table:
  `atom_index, id, x, y, z, bucket, structure_id (0/1/2), structure_name, cluster_id`
  plus per-l scalar columns `q_<l>` (e.g. `q_4 … q_12`) and the stable `q6`
  classifier column, the full q_l set as a single `q` `list<double>` column, and
  (when enabled) `w_<l>`, `w_hat_<l>`, and a `components` `list<double>` column.
  `bucket`/`structure_name` is the q_6-based heuristic label
  (`Liquid-like`, `Interface`, `Crystal-like`) driving the AtomisticExporter GLB.

