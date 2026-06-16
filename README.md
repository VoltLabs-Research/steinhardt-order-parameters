# Steinhardt Order Parameters

Computes per-atom bond-orientational order parameters q_l, Wigner-3j invariants w_l and normalized w_l_hat (ports LAMMPS `compute orientorder/atom`) for crystallinity, nucleation and liquid/solid discrimination.

## Install

```bash
vpm install @voltlabs/steinhardt-order-parameters
```

## CLI

```bash
steinhardt-order-parameters <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--qlist <csv>` | no | `4,6,8,10,12` | Comma-separated list of orders l to compute. |
| `--nnn <int>` | no | `12` | Fixed number of nearest neighbors (`0` = use `--cutoff`). |
| `--cutoff <float>` | no | `0.0` | Cutoff radius (only used when `--nnn 0`). |
| `--wl <bool>` | no | `false` | Also compute w_l (Wigner-3j invariant). |
| `--wlHat <bool>` | no | `false` | Also compute w_l_hat = w_l / q_l^3 (dimensionless). |
| `--components <l>` | no | `-1` | Emit the 2l+1 complex Q_l^m/\|Q_l\| components for a single l (disabled if < 0). |
| `--onlySelected` | no | `false` | Only classify atoms flagged by a `Selection` column. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_steinhardt.parquet` | Steinhardt Order | — (listing-only) |
| `{output_base}_atoms.parquet` | Steinhardt Model | AtomisticExporter → glb |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins
