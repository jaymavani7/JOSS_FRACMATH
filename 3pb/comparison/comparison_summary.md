# Abaqus vs MATLAB 3PB Comparison Summary

This report compares the computational results and performance of the three-point bending (3PB) simulation using Abaqus (UMAT) and the MATLAB solver (FRACMATH).

## Summary Table (Markdown)

| Quantity | MATLAB | Abaqus + UMAT | Ratio |
| :--- | :---: | :---: | :---: |
| **Peak load (kN)** | 3.90 | 3.84 | 1.017 |
| **CMOD at peak (mm)** | 0.026 | 0.026 | 0.978 |
| **Wall-clock time (s)** | 97.18 | 1336.00 | 0.073 |
| **Total process time (s)** | 103.90 | 1356.20 | 0.077 |

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
Peak load (kN) & 3.90 & 3.84 & 1.017 \\
CMOD at peak (mm) & 0.026 & 0.026 & 0.978 \\
Wall-clock time (s) & 97.18 & 1336.00 & 0.073 \\
Total process (s) & 103.90 & 1356.20 & 0.077 \\
\hline
\end{tabular}
\end{table}
```

## Visualizations
- Combined bar and breakdown pie chart: [performance_breakdown.png](performance_breakdown.png)
- Load vs CMOD comparison plot: [comparison_load_cmod.png](comparison_load_cmod.png)
