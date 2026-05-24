import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

def generate_comparison():
    # Paths
    matlab_csv = r'Matlab\Gregoire_3PB\results\matlab_load_cmod.csv'
    abaqus_csv = r'abaqus\results\load_cmod_data.csv'
    matlab_timing_file = r'Matlab\Gregoire_3PB\results\matlab_timing.txt'
    
    # Output dirs
    plot_dir = r'Comparison\plot'
    time_dir = r'Comparison\time'
    mem_dir = r'Comparison\memory'
    
    # Load data
    m_data = pd.read_csv(matlab_csv)
    a_data = pd.read_csv(abaqus_csv)
    
    # Extract data for plotting
    # Matlab: CMOD [mm], Load [N]
    # Abaqus: time, load_N, cmod_mm
    
    m_cmod = m_data.iloc[:, 0]
    m_load = m_data.iloc[:, 1] / 1000.0 # Convert to kN
    
    a_cmod = a_data['cmod_mm']
    a_load = a_data['load_N'] / 1000.0 # Convert to kN
    
    # Calculate peaks
    m_peak_load = m_load.max()
    m_cmod_at_peak = m_cmod[m_load.idxmax()]
    
    a_peak_load = a_load.max()
    a_cmod_at_peak = a_cmod[a_load.idxmax()]
    
    # --- Plotting ---
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": '--'
    })
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    
    # Plot Matlab
    ax.plot(m_cmod, m_load, label='Matlab (Reference)', color='#1f77b4', linewidth=2, zorder=2)
    
    # Plot Abaqus (with markers)
    ax.plot(a_cmod, a_load, label='Abaqus + UMAT', color='#ff7f0e', linestyle='--', linewidth=1.5, marker='o', markevery=10, markerfacecolor='white', markeredgewidth=1.5, zorder=3)
    
    ax.set_xlabel('CMOD [mm]', fontweight='bold')
    ax.set_ylabel('Load [kN]', fontweight='bold')
    ax.set_title('Load vs. CMOD: Matlab vs. Abaqus Comparison', fontweight='bold', pad=15)
    ax.legend(frameon=True, fancybox=True, shadow=True)
    
    # Set limits with some padding
    ax.set_xlim(0, max(m_cmod.max(), a_cmod.max()) * 1.05)
    ax.set_ylim(0, max(m_load.max(), a_load.max()) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'load_cmod_comparison.png'), dpi=300)
    plt.savefig(os.path.join(plot_dir, 'load_cmod_comparison.pdf'))
    plt.close()
    
    # --- Bar Charts for Timing and Memory ---
    # Data
    labels = ['Matlab', 'Abaqus']
    times = [5.00, 19.00]
    mems = [30.9, 32.0]
    
    # Timing Bar Chart
    fig, ax = plt.subplots(figsize=(6, 5), dpi=150)
    bars = ax.bar(labels, times, color=['#1f77b4', '#ff7f0e'], alpha=0.8, edgecolor='black', linewidth=1.2)
    ax.set_ylabel('Wall-clock Time [s]', fontweight='bold')
    ax.set_title('Computational Time Comparison', fontweight='bold', pad=15)
    
    # Add values on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}s',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'time_comparison_bar.png'), dpi=300)
    plt.close()
    
    # Memory Bar Chart
    fig, ax = plt.subplots(figsize=(6, 5), dpi=150)
    bars = ax.bar(labels, mems, color=['#1f77b4', '#ff7f0e'], alpha=0.8, edgecolor='black', linewidth=1.2)
    ax.set_ylabel('Peak RAM [MB]', fontweight='bold')
    ax.set_title('Memory Usage Comparison', fontweight='bold', pad=15)
    ax.set_ylim(0, max(mems) * 1.2)
    
    # Add values on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.1f} MB',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'memory_comparison_bar.png'), dpi=300)
    plt.close()
    
    # --- Timing and Memory Reports ---
    timing_report = f"""MATLAB VS ABAQUS TIMING COMPARISON
========================================
Quantity            Matlab      Abaqus      Ratio (A/M)
----------------------------------------
Wall-clock (s)      5.00        19.00       {19.00/5.00:.2f}x
Load Steps          200         200         1.00x
DOFs                5808        5808        1.00x
Peak Load (kN)      {m_peak_load:.3f}       {a_peak_load:.3f}       {a_peak_load/m_peak_load:.4f}x
CMOD@Peak (mm)      {m_cmod_at_peak:.5f}     {a_cmod_at_peak:.5f}     {a_cmod_at_peak/m_cmod_at_peak:.4f}x
----------------------------------------
"""
    
    memory_report = f"""MATLAB VS ABAQUS MEMORY COMPARISON
========================================
Quantity            Matlab      Abaqus      Ratio (A/M)
----------------------------------------
Peak RAM (MB)       30.9        32.0        {32.0/30.9:.2f}x
----------------------------------------
Note: Abaqus memory is based on the 'Memory to minimize I/O' estimate.
"""
    
    with open(os.path.join(time_dir, 'timing_comparison.txt'), 'w') as f:
        f.write(timing_report)
        
    with open(os.path.join(mem_dir, 'memory_comparison.txt'), 'w') as f:
        f.write(memory_report)
        
    print("Comparison generated successfully.")

if __name__ == "__main__":
    generate_comparison()
