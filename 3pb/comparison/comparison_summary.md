# Abaqus vs MATLAB 3PB Comparison Summary

This report compares the computational results and performance of the three-point bending (3PB) simulation using Abaqus (UMAT) and the MATLAB solver (FRACMATH).

## Summary Table (Markdown)

| Quantity | MATLAB | Abaqus + UMAT | Ratio |
| :--- | :---: | :---: | :---: |
| **Peak load (kN)** | 3.64 | 3.61 | 1.009 |
| **CMOD at peak (mm)** | 0.023 | 0.022 | 1.014 |
| **Solver/submit wall-clock (s)** | 547.58 | 1996.25 | 0.274 |
| **End-to-end / internal time (s)** | 590.54 | 1968.00 | 0.300 |
| **Load-CMOD rows** | 10000 | 10001 | - |

## LaTeX Table Format
```latex
\begin{table}[htbp]
\centering
\caption{MATLAB vs. Abaqus on the 2D 3PB case (identical mesh, material, and tolerances). Abaqus ran on 4 threads, MATLAB on 1.}
\label{tab:comparison}
\begin{tabular}{lrrr}
\hline
Quantity & MATLAB & Abaqus + UMAT & Ratio \\
\hline
Peak load (kN) & 3.64 & 3.61 & 1.009 \\
CMOD at peak (mm) & 0.023 & 0.022 & 1.014 \\
Solver/submit wall-clock (s) & 547.58 & 1996.25 & 0.274 \\
End-to-end / internal time (s) & 590.54 & 1968.00 & 0.300 \\
\hline
\end{tabular}
\end{table}
```

## Visualizations
- Combined bar and breakdown pie chart: [performance_breakdown.png](performance_breakdown.png)
- Load vs CMOD comparison plot: [comparison_load_cmod.png](comparison_load_cmod.png)
- Side-by-side damage/mesh panels: [comparison_damage_last_step.png](comparison_damage_last_step.png)
