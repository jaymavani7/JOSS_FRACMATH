# Theory manual

This folder contains the theory manual for the FRACMATH formulation and solver implementation.

## Files

| File | Purpose |
| --- | --- |
| `theory_manual.tex` | LaTeX source |
| `theory_manual.pdf` | Compiled PDF manual |

## What the manual covers

- Scalar isotropic continuum damage mechanics.
- Modified von Mises equivalent strain.
- Exponential softening.
- Crack-band regularization.
- Direction-dependent Oliver bandwidth.
- Newton-Raphson solution procedure.
- Notes connecting the MATLAB and Abaqus/UMAT implementations.

## Rebuild command

From this folder, run:

```bash
pdflatex theory_manual.tex
pdflatex theory_manual.tex
```

A LaTeX distribution is required only if you want to rebuild the PDF.
