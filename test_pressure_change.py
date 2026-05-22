import matplotlib.pyplot as plt
import numpy as np
from previous_HXs import hxs, mass_flows, pressure_changes
from hydraulic_analysis import solve_mass_flows

# -----------------------------
# Collect data
# -----------------------------
tube_lengths = []
tube_velocities_hot = []
tube_velocities_cold = []
num_tubes = []
num_baffles = []
num_tube_passes = []
num_shell_passes = []
rel_shell_passses = []
hot_errors = []
cold_errors = []
rel_cold_errors = []
rel_hot_errors = []
dp_hot_measured_list = []
dp_cold_measured_list = []
dp_hot_pred_list = []
dp_cold_pred_list = []

for i, hx in enumerate(hxs):
#     m_cold, m_hot = mass_flows[i]
    m_cold, m_hot = solve_mass_flows(hx)

    # Predictions
    dp_cold_pred = hx.cold_side_pressure_drop(m_cold)
    dp_hot_pred = hx.hot_side_pressure_drop(m_hot)
    
    # Experimental values
    dp_cold_measured, dp_hot_measured = pressure_changes[i]
    
    # Velocities
    v_hot = hx.tube_side_velocity(m_hot)
    v_cold = hx.shell_side_velocity(m_cold)
    
    # Absolute errors (Pa)
    abs_error_cold = dp_cold_pred - dp_cold_measured
    abs_error_hot = dp_hot_pred - dp_hot_measured

    rel_error_cold = abs_error_cold / dp_cold_measured
    rel_cold_errors.append(rel_error_cold)

    rel_error_hot = abs_error_cold / dp_hot_measured
    rel_hot_errors.append(rel_error_hot)

    
    # Store parameters
    tube_lengths.append(hx.tube_length)
    tube_velocities_hot.append(v_hot)
    tube_velocities_cold.append(v_cold)
    num_tubes.append(hx.number_of_tubes)
    num_baffles.append(hx.number_of_baffles)
    num_tube_passes.append(hx.tube_passes)
    num_shell_passes.append(hx.shell_passes)
    rel_shell_passses.append(hx.tube_passes/hx.shell_passes)
    
    # Store errors and pressures
    hot_errors.append(abs_error_hot)
    cold_errors.append(abs_error_cold)
    dp_hot_measured_list.append(dp_hot_measured)
    dp_cold_measured_list.append(dp_cold_measured)
    dp_hot_pred_list.append(dp_hot_pred)
    dp_cold_pred_list.append(dp_cold_pred)

#     print(hx.cold_side_contributions)
#     print(hx.shell_side_reynolds_number(m_cold))

# Convert to numpy arrays for easier manipulation
tube_lengths = np.array(tube_lengths)
tube_velocities_hot = np.array(tube_velocities_hot)
tube_velocities_cold = np.array(tube_velocities_cold)
num_tubes = np.array(num_tubes)
num_baffles = np.array(num_baffles)
num_tube_passes = np.array(num_tube_passes)
hot_errors = np.array(hot_errors)
cold_errors = np.array(cold_errors)

# -----------------------------
# HOT SIDE PLOTS
# -----------------------------
# fig, axes = plt.subplots(2, 2, figsize=(14, 10))
# fig.suptitle('HOT SIDE (Tube) Pressure Drop Error Analysis', fontsize=16, fontweight='bold')

# # Plot 1: Error vs Tube Length
# ax = axes[0, 0]
# ax.scatter(tube_lengths, hot_errors, s=100, alpha=0.6, edgecolors='k', c='orange')
# ax.axhline(0, color='r', linestyle='--', alpha=0.5, label='Zero error')
# # Linear fit
# m, c = np.polyfit(tube_lengths, hot_errors, 1)
# xfit = np.linspace(min(tube_lengths), max(tube_lengths), 100)
# ax.plot(xfit, m*xfit + c, 'b-', linewidth=2, label=f'Fit: y={m:.1f}x+{c:.1f}')
# ax.set_xlabel('Tube Length (m)', fontsize=12)
# ax.set_ylabel('Pressure Error (Pa)', fontsize=12)
# ax.set_title('Error vs Tube Length', fontsize=13, fontweight='bold')
# ax.grid(True, alpha=0.3)
# ax.legend()
# # Add correlation coefficient
# r = np.corrcoef(tube_lengths, hot_errors)[0, 1]
# ax.text(0.05, 0.95, f'R = {r:.3f}', transform=ax.transAxes, 
#         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# # Plot 2: Error vs Tube Velocity
# ax = axes[0, 1]
# ax.scatter(tube_velocities_hot, hot_errors, s=100, alpha=0.6, edgecolors='k', c='orange')
# ax.axhline(0, color='r', linestyle='--', alpha=0.5, label='Zero error')
# # Linear fit
# m, c = np.polyfit(tube_velocities_hot, hot_errors, 1)
# xfit = np.linspace(min(tube_velocities_hot), max(tube_velocities_hot), 100)
# ax.plot(xfit, m*xfit + c, 'b-', linewidth=2, label=f'Fit: y={m:.1f}x+{c:.1f}')
# ax.set_xlabel('Tube Velocity (m/s)', fontsize=12)
# ax.set_ylabel('Pressure Error (Pa)', fontsize=12)
# ax.set_title('Error vs Tube Velocity', fontsize=13, fontweight='bold')
# ax.grid(True, alpha=0.3)
# ax.legend()
# r = np.corrcoef(tube_velocities_hot, hot_errors)[0, 1]
# ax.text(0.05, 0.95, f'R = {r:.3f}', transform=ax.transAxes, 
#         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# # Plot 3: Error vs Number of Tubes
# ax = axes[1, 0]
# ax.scatter(num_tubes, hot_errors, s=100, alpha=0.6, edgecolors='k', c='orange')
# ax.axhline(0, color='r', linestyle='--', alpha=0.5, label='Zero error')
# # Linear fit
# m, c = np.polyfit(num_tubes, hot_errors, 1)
# xfit = np.linspace(min(num_tubes), max(num_tubes), 100)
# ax.plot(xfit, m*xfit + c, 'b-', linewidth=2, label=f'Fit: y={m:.1f}x+{c:.1f}')
# ax.set_xlabel('Number of Tubes', fontsize=12)
# ax.set_ylabel('Pressure Error (Pa)', fontsize=12)
# ax.set_title('Error vs Number of Tubes', fontsize=13, fontweight='bold')
# ax.grid(True, alpha=0.3)
# ax.legend()
# r = np.corrcoef(num_tubes, hot_errors)[0, 1]
# ax.text(0.05, 0.95, f'R = {r:.3f}', transform=ax.transAxes, 
#         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# # Plot 4: Error vs Number of Tube Passes
# ax = axes[1, 1]
# # Get unique tube passes for better visualization
# unique_passes = np.unique(num_tube_passes)
# for tp in unique_passes:
#     mask = num_tube_passes == tp
#     ax.scatter(num_tube_passes[mask], hot_errors[mask], s=100, alpha=0.6, 
#               edgecolors='k', label=f'{tp} passes')
# ax.axhline(0, color='r', linestyle='--', alpha=0.5, label='Zero error')
# ax.set_xlabel('Number of Tube Passes', fontsize=12)
# ax.set_ylabel('Pressure Error (Pa)', fontsize=12)
# ax.set_title('Error vs Number of Tube Passes', fontsize=13, fontweight='bold')
# ax.grid(True, alpha=0.3)
# ax.legend()
# if len(unique_passes) > 1:
#     r = np.corrcoef(num_tube_passes, hot_errors)[0, 1]
#     ax.text(0.05, 0.95, f'R = {r:.3f}', transform=ax.transAxes, 
#             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# plt.tight_layout()
# plt.savefig('hot_side_pressure_error_analysis.png', dpi=300, bbox_inches='tight')

# -----------------------------
# COLD SIDE PLOTS
# -----------------------------

errors = [a for a in rel_hot_errors]
hotcold = 'HOT'
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(hotcold + ' SIDE (Shell) Pressure Drop Error Analysis', fontsize=16, fontweight='bold')

# Plot 1: Error vs Shell Passes
x1 = tube_lengths
label1 = 'Tube Lengths'

ax = axes[0, 0]
ax.scatter(x1, errors, s=100, alpha=0.6, edgecolors='k', c='cyan')
ax.axhline(0, color='r', linestyle='--', alpha=0.5, label='Zero error')
# Linear fit
m, c = np.polyfit(x1, errors, 1)
xfit = np.linspace(min(x1), max(x1), 100)
# ax.plot(xfit, m*xfit + c, 'b-', linewidth=2, label=f'Fit: y={m:.1f}x+{c:.1f}')
ax.set_xlabel(f'{label1}', fontsize=12)
ax.set_ylabel('Pressure Error', fontsize=12)
ax.set_title(f'Error vs {label1}', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()
r = np.corrcoef(x1, errors)[0, 1]
ax.text(0.05, 0.95, f'R = {r:.3f}', transform=ax.transAxes, 
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Plot 2: Error vs Shell Velocity
x2 = tube_velocities_hot
label2 = 'Velocity Hot (m/s)'

ax = axes[0, 1]
ax.scatter(x2, errors, s=100, alpha=0.6, edgecolors='k', c='cyan')
ax.axhline(0, color='r', linestyle='--', alpha=0.5, label='Zero error')
# Linear fit
m, c = np.polyfit(tube_velocities_cold, errors, 1)
xfit = np.linspace(min(tube_velocities_cold), max(tube_velocities_cold), 100)
# ax.plot(xfit, m*xfit + c, 'b-', linewidth=2, label=f'Fit: y={m:.1f}x+{c:.1f}')
ax.set_xlabel(label2, fontsize=12)
ax.set_ylabel('Pressure Error', fontsize=12)
ax.set_title('Error vs Shell Velocity', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()
r = np.corrcoef(tube_velocities_cold, errors)[0, 1]
ax.text(0.05, 0.95, f'R = {r:.3f}', transform=ax.transAxes, 
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Plot 3: Error vs Number of Tubes
ax = axes[1, 0]
ax.scatter(num_tubes, errors, s=100, alpha=0.6, edgecolors='k', c='cyan')
ax.axhline(0, color='r', linestyle='--', alpha=0.5, label='Zero error')
# Linear fit
m, c = np.polyfit(num_tubes, errors, 1)
xfit = np.linspace(min(num_tubes), max(num_tubes), 100)
# ax.plot(xfit, m*xfit + c, 'b-', linewidth=2, label=f'Fit: y={m:.1f}x+{c:.1f}')
ax.set_xlabel('Number of Tubes', fontsize=12)
ax.set_ylabel('Pressure Error', fontsize=12)
ax.set_title('Error vs Number of Tubes', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()
r = np.corrcoef(num_tubes, errors)[0, 1]
ax.text(0.05, 0.95, f'R = {r:.3f}', transform=ax.transAxes, 
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Plot 4: Error vs Number of Baffles
ax = axes[1, 1]
ax.scatter(num_baffles, errors, s=100, alpha=0.6, edgecolors='k', c='cyan')
ax.axhline(0, color='r', linestyle='--', alpha=0.5, label='Zero error')
# Linear fit
m, c = np.polyfit(num_baffles, errors, 1)
xfit = np.linspace(min(num_baffles), max(num_baffles), 100)
# ax.plot(xfit, m*xfit + c, 'b-', linewidth=2, label=f'Fit: y={m:.1f}x+{c:.1f}')
ax.set_xlabel('Number of Baffles', fontsize=12)
ax.set_ylabel('Pressure Error', fontsize=12)
ax.set_title('Error vs Number of Baffles', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()
r = np.corrcoef(num_baffles, errors)[0, 1]
ax.text(0.05, 0.95, f'R = {r:.3f}', transform=ax.transAxes, 
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('cold_side_pressure_error_analysis.png', dpi=300, bbox_inches='tight')

# -----------------------------
# SUMMARY STATISTICS
# -----------------------------
print("\n" + "="*60)
print("PRESSURE DROP ERROR ANALYSIS")
print("="*60)

# print("\nHOT SIDE (Tube) Statistics:")
# print(f"  Mean error:              {np.mean(hot_errors):+.1f} Pa")
# print(f"  Std dev:                 {np.std(hot_errors):.1f} Pa")
# print(f"  Max overprediction:      {np.max(hot_errors):+.1f} Pa")
# print(f"  Max underprediction:     {np.min(hot_errors):+.1f} Pa")
# print(f"  RMSE:                    {np.sqrt(np.mean(hot_errors**2)):.1f} Pa")
# print(f"  Mean absolute error:     {np.mean(np.abs(hot_errors)):.1f} Pa")

print("\nCOLD SIDE (Shell) Statistics:")
print(f"  Mean error:              {np.mean(cold_errors):+.1f} Pa")
print(f"  Std dev:                 {np.std(cold_errors):.1f} Pa")
print(f"  Max overprediction:      {np.max(cold_errors):+.1f} Pa")
print(f"  Max underprediction:     {np.min(cold_errors):+.1f} Pa")
print(f"  RMSE:                    {np.sqrt(np.mean(cold_errors**2)):.1f} Pa")
print(f"  Mean absolute error:     {np.mean(np.abs(cold_errors)):.1f} Pa")

# print("\nCorrelation Analysis (HOT SIDE):")
# print(f"  Error vs Tube Length:    R = {np.corrcoef(tube_lengths, hot_errors)[0,1]:.3f}")
# print(f"  Error vs Tube Velocity:  R = {np.corrcoef(tube_velocities_hot, hot_errors)[0,1]:.3f}")
# print(f"  Error vs Num Tubes:      R = {np.corrcoef(num_tubes, hot_errors)[0,1]:.3f}")
# print(f"  Error vs Tube Passes:    R = {np.corrcoef(num_tube_passes, hot_errors)[0,1]:.3f}")

print("\nCorrelation Analysis (COLD SIDE):")
print(f"  Error vs Tube Length:    R = {np.corrcoef(tube_lengths, cold_errors)[0,1]:.3f}")
print(f"  Error vs Shell Velocity: R = {np.corrcoef(tube_velocities_cold, cold_errors)[0,1]:.3f}")
print(f"  Error vs Num Tubes:      R = {np.corrcoef(num_tubes, cold_errors)[0,1]:.3f}")
print(f"  Error vs Num Baffles:    R = {np.corrcoef(num_baffles, cold_errors)[0,1]:.3f}")

plt.show()