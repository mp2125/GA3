import numpy as np
from scipy.optimize import minimize, differential_evolution
import matplotlib.pyplot as plt
from copy import deepcopy

# You'll need to import your actual modules
from previous_HXs import hxs, mass_flows, pressure_changes


class HeatExchangerOptimizer:
    """
    Calibration and diagnostic tool for shell-and-tube heat exchanger model
    """
    
    def __init__(self, hxs, mass_flows, pressure_changes):
        """
        Args:
            hxs: List of ShellAndTubeHeatExchanger instances
            mass_flows: List of (m_dot_cold, m_dot_hot) tuples
            pressure_changes: List of (dp_cold_measured, dp_hot_measured) tuples
        """
        self.hxs = hxs
        self.mass_flows = mass_flows
        self.pressure_changes = pressure_changes
        self.original_hxs = [deepcopy(hx) for hx in hxs]
        
    def diagnose_single_hx(self, idx, verbose=True):
        """
        Detailed breakdown of pressure drop components for one heat exchanger
        """
        hx = self.hxs[idx]
        m_cold, m_hot = self.mass_flows[idx]
        dp_cold_meas, dp_hot_meas = self.pressure_changes[idx]
        
        # ===== HOT SIDE (TUBE) BREAKDOWN =====
        velocity_hot = hx.tube_side_velocity(m_hot)
        reynolds_hot = hx.tube_side_reynolds_number(velocity_hot)
        
        relative_roughness = hx.tube_roughness / hx.tube_inner_diameter
        friction_factor = hx.friction_factor(reynolds_hot, relative_roughness)
        friction_factor_used = friction_factor / hx.friction_divisor
        
        effective_length = hx.tube_length * hx.tube_passes
        pipe_friction = (friction_factor_used * (effective_length / hx.tube_inner_diameter) 
                        * 0.5 * hx.fluid_density * velocity_hot**2)
        
        Kc, Ke = hx.get_entrance_exit_coefficients(reynolds_hot)
        entrance_exit = (hx.tube_passes * 0.5 * hx.fluid_density 
                         * velocity_hot**2 * (Kc + Ke))
        
        nozzle_velocity_hot = m_hot / (hx.fluid_density * hx.nozzle_area_tube_side)
        nozzle_hot = 2 * hx.nozzle_correction_factor * 0.5 * hx.fluid_density * nozzle_velocity_hot**2
        
        hose_hot = 2 * hx.hose_pressure_drop(m_hot, hx.K_hose_hot)
        
        misc_hot = hx.K_tube_misc * 0.5 * hx.fluid_density * velocity_hot**2
        
        total_hot = pipe_friction + entrance_exit + nozzle_hot + hose_hot + misc_hot
        
        # ===== COLD SIDE (SHELL) BREAKDOWN =====
        velocity_cold = hx.shell_side_velocity(m_cold)
        reynolds_cold = hx.shell_side_reynolds_number(velocity_cold)
        a = hx.shell_side_friction_coefficient(reynolds_cold)
        N = hx.number_of_baffles + 1
        
        bundle_dp = 4 * a * N * hx.fluid_density * velocity_cold**2 * hx.shell_passes
        
        # Turning losses at baffles
        turning_dp = (hx.number_of_baffles * hx.K_turn * 0.5 * 
                      hx.fluid_density * velocity_cold**2 * hx.shell_passes)
        
        nozzle_velocity_cold = m_cold / (hx.fluid_density * hx.nozzle_area_shell_side)
        nozzle_cold = 2 * hx.nozzle_correction_factor * 0.5 * hx.fluid_density * nozzle_velocity_cold**2
        
        hose_cold = 2 * hx.hose_pressure_drop(m_cold, hx.K_hose_cold)
        
        total_cold = bundle_dp + turning_dp + nozzle_cold + hose_cold
        
        if verbose:
            print(f"\n{'='*60}")
            print(f"HEAT EXCHANGER #{idx + 1}")
            print(f"{'='*60}")
            print(f"Geometry: {hx.number_of_tubes} tubes, {hx.number_of_baffles} baffles, L={hx.tube_length:.3f}m")
            print(f"Mass flows: cold={m_cold:.3f} kg/s, hot={m_hot:.3f} kg/s")
            
            print(f"\n--- HOT SIDE (TUBE) ---")
            print(f"  Velocity:           {velocity_hot:.3f} m/s")
            print(f"  Reynolds:           {reynolds_hot:.0f}")
            print(f"  Friction factor:    {friction_factor:.6f} (using {friction_factor_used:.6f})")
            print(f"  Friction divisor:   {hx.friction_divisor:.2f}")
            print(f"  Kc, Ke:             {Kc:.4f}, {Ke:.4f}")
            print(f"  Entrance/exit mult: {hx.entrance_exit_multiplier:.2f}")
            print(f"\n  Pipe friction:      {pipe_friction:8.1f} Pa  ({100*pipe_friction/total_hot:5.1f}%)")
            print(f"  Entrance/exit:      {entrance_exit:8.1f} Pa  ({100*entrance_exit/total_hot:5.1f}%)")
            print(f"  Nozzles:            {nozzle_hot:8.1f} Pa  ({100*nozzle_hot/total_hot:5.1f}%)")
            print(f"  Hoses:              {hose_hot:8.1f} Pa  ({100*hose_hot/total_hot:5.1f}%)")
            print(f"  Misc (K={hx.K_tube_misc:.2f}):      {misc_hot:8.1f} Pa  ({100*misc_hot/total_hot:5.1f}%)")
            print(f"  ---")
            print(f"  TOTAL PREDICTED:    {total_hot:8.1f} Pa")
            print(f"  MEASURED:           {dp_hot_meas:8.1f} Pa")
            print(f"  ERROR:              {100*(total_hot-dp_hot_meas)/dp_hot_meas:6.1f}%")
            
            print(f"\n--- COLD SIDE (SHELL) ---")
            print(f"  Velocity:           {velocity_cold:.3f} m/s")
            print(f"  Reynolds:           {reynolds_cold:.0f}")
            print(f"  Friction coeff (a): {a:.4f}")
            print(f"  Shell friction mult:{hx.shell_friction_a / (0.34 if hx.is_square_layout else 0.2):.2f}")
            print(f"  Crossflow corr:     {hx.crossflow_correction_factor:.2f}")
            print(f"  K_turn:             {hx.K_turn:.2f}")
            print(f"\n  Bundle crossflow:   {bundle_dp:8.1f} Pa  ({100*bundle_dp/total_cold:5.1f}%)")
            print(f"  Turning losses:     {turning_dp:8.1f} Pa  ({100*turning_dp/total_cold:5.1f}%)")
            print(f"  Nozzles:            {nozzle_cold:8.1f} Pa  ({100*nozzle_cold/total_cold:5.1f}%)")
            print(f"  Hoses:              {hose_cold:8.1f} Pa  ({100*hose_cold/total_cold:5.1f}%)")
            print(f"  ---")
            print(f"  TOTAL PREDICTED:    {total_cold:8.1f} Pa")
            print(f"  MEASURED:           {dp_cold_meas:8.1f} Pa")
            print(f"  ERROR:              {100*(total_cold-dp_cold_meas)/dp_cold_meas:6.1f}%")
        
        return {
            'hot': {
                'pipe_friction': pipe_friction,
                'entrance_exit': entrance_exit,
                'nozzle': nozzle_hot,
                'hose': hose_hot,
                'misc': misc_hot,
                'total_predicted': total_hot,
                'measured': dp_hot_meas,
                'error_pct': 100*(total_hot-dp_hot_meas)/dp_hot_meas,
                'reynolds': reynolds_hot,
                'friction_factor': friction_factor
            },
            'cold': {
                'bundle': bundle_dp,
                'turning': turning_dp,
                'nozzle': nozzle_cold,
                'hose': hose_cold,
                'total_predicted': total_cold,
                'measured': dp_cold_meas,
                'error_pct': 100*(total_cold-dp_cold_meas)/dp_cold_meas,
                'reynolds': reynolds_cold
            }
        }
    
    def diagnose_all(self):
        """Run diagnostics on all heat exchangers"""
        results = []
        for i in range(len(self.hxs)):
            results.append(self.diagnose_single_hx(i, verbose=True))
        return results
    
    def objective_function(self, params, param_names, weight_cold=1.0, weight_hot=1.0):
        """
        Objective function for optimization
        
        Args:
            params: Array of parameter values to optimize
            param_names: List of parameter names being optimized
            weight_cold: Weight for cold-side errors
            weight_hot: Weight for hot-side errors
        
        Returns:
            Total weighted sum of squared relative errors
        """
        # Apply parameters to all heat exchangers
        for hx in self.hxs:
            for param_name, param_value in zip(param_names, params):
                if param_name == 'K_tube_misc':
                    hx.K_tube_misc = param_value
                elif param_name == 'crossflow_correction_factor':
                    hx.crossflow_correction_factor = param_value
                elif param_name == 'shell_friction_multiplier':
                    # Multiply the base value
                    base_a = 0.34 if hx.is_square_layout else 0.2
                    hx.shell_friction_a = base_a * param_value
                elif param_name == 'nozzle_correction_factor':
                    hx.nozzle_correction_factor = param_value
                elif param_name == 'friction_divisor':
                    hx.friction_divisor = param_value
                elif param_name == 'entrance_exit_multiplier':
                    hx.entrance_exit_multiplier = param_value
                elif param_name == 'K_turn':
                    hx.K_turn = param_value
        
        # Calculate total error
        total_error = 0
        for i, hx in enumerate(self.hxs):
            m_cold, m_hot = self.mass_flows[i]
            dp_cold_meas, dp_hot_meas = self.pressure_changes[i]
            
            # Get predictions
            dp_cold_pred = hx.cold_side_pressure_drop(m_cold)
            dp_hot_pred = hx.hot_side_pressure_drop(m_hot)
            
            # Relative errors
            error_cold = ((dp_cold_pred - dp_cold_meas) / dp_cold_meas) ** 2
            error_hot = ((dp_hot_pred - dp_hot_meas) / dp_hot_meas) ** 2
            
            total_error += weight_cold * error_cold + weight_hot * error_hot
        
        return total_error
    
    def optimize_parameters(self, param_config, method='differential_evolution', 
                          weight_cold=1.0, weight_hot=1.0):
        """
        Optimize correction factors
        
        Args:
            param_config: Dict with parameter names as keys and (min, max, initial) as values
                Example: {
                    'K_tube_misc': (0.0, 5.0, 0.0),
                    'crossflow_correction_factor': (0.5, 5.0, 1.6),
                    'shell_friction_multiplier': (0.5, 3.0, 1.3),
                    'nozzle_correction_factor': (0.1, 2.0, 0.5),
                    'friction_divisor': (0.5, 10.0, 1.0),
                    'entrance_exit_multiplier': (0.1, 2.0, 0.5),
                    'K_turn': (0.0, 5.0, 1.0)
                }
            method: 'differential_evolution' (global) or 'minimize' (local)
            weight_cold: Weight for cold-side errors
            weight_hot: Weight for hot-side errors
        
        Returns:
            Optimization result object
        """
        param_names = list(param_config.keys())
        bounds = [param_config[name][:2] for name in param_names]
        x0 = [param_config[name][2] for name in param_names]
        
        print(f"\nOptimizing {len(param_names)} parameters:")
        print(f"Weights: cold={weight_cold:.2f}, hot={weight_hot:.2f}")
        for name, (lb, ub, init) in param_config.items():
            print(f"  {name:30s}: [{lb:.3f}, {ub:.3f}], initial={init:.3f}")
        
        if method == 'differential_evolution':
            result = differential_evolution(
                lambda x: self.objective_function(x, param_names, weight_cold, weight_hot),
                bounds=bounds,
                maxiter=200,
                popsize=15,
                seed=42,
                disp=True,
                atol=1e-8,
                tol=1e-8,
                workers=1
            )
        else:
            result = minimize(
                lambda x: self.objective_function(x, param_names, weight_cold, weight_hot),
                x0=x0,
                bounds=bounds,
                method='L-BFGS-B',
                options={'maxiter': 1000}
            )
        
        print(f"\n{'='*60}")
        print("OPTIMIZATION RESULTS")
        print(f"{'='*60}")
        print(f"Final objective value: {result.fun:.6f}")
        print(f"\nOptimal parameters:")
        for name, value in zip(param_names, result.x):
            print(f"  {name:30s} = {value:.6f}")
        
        # Apply optimal parameters
        self.objective_function(result.x, param_names, weight_cold, weight_hot)
        
        return result
    
    def plot_comparison(self, save_path=None):
        """
        Create comparison plots of predicted vs measured pressure drops
        """
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        cold_predicted = []
        cold_measured = []
        hot_predicted = []
        hot_measured = []
        
        for i, hx in enumerate(self.hxs):
            m_cold, m_hot = self.mass_flows[i]
            dp_cold_meas, dp_hot_meas = self.pressure_changes[i]
            
            dp_cold_pred = hx.cold_side_pressure_drop(m_cold)
            dp_hot_pred = hx.hot_side_pressure_drop(m_hot)
            
            cold_predicted.append(dp_cold_pred)
            cold_measured.append(dp_cold_meas)
            hot_predicted.append(dp_hot_pred)
            hot_measured.append(dp_hot_meas)
        
        # Cold side
        ax = axes[0]
        ax.scatter(cold_measured, cold_predicted, s=100, alpha=0.6, edgecolors='k')
        max_val = max(max(cold_measured), max(cold_predicted))
        min_val = min(min(cold_measured), min(cold_predicted))
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='Perfect prediction')
        ax.plot([min_val, max_val], [0.9*min_val, 0.9*max_val], 'g:', alpha=0.5, label='±10%')
        ax.plot([min_val, max_val], [1.1*min_val, 1.1*max_val], 'g:', alpha=0.5)
        ax.set_xlabel('Measured ΔP (Pa)', fontsize=12)
        ax.set_ylabel('Predicted ΔP (Pa)', fontsize=12)
        ax.set_title('Cold Side (Shell)', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Calculate R² and RMSE
        cold_r2 = 1 - np.sum((np.array(cold_measured) - np.array(cold_predicted))**2) / \
                     np.sum((np.array(cold_measured) - np.mean(cold_measured))**2)
        cold_rmse = np.sqrt(np.mean((np.array(cold_measured) - np.array(cold_predicted))**2))
        cold_mape = np.mean(np.abs((np.array(cold_measured) - np.array(cold_predicted)) / np.array(cold_measured))) * 100
        
        ax.text(0.05, 0.95, f'R² = {cold_r2:.4f}\nRMSE = {cold_rmse:.1f} Pa\nMAPE = {cold_mape:.1f}%',
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Hot side
        ax = axes[1]
        ax.scatter(hot_measured, hot_predicted, s=100, alpha=0.6, edgecolors='k', color='orange')
        max_val = max(max(hot_measured), max(hot_predicted))
        min_val = min(min(hot_measured), min(hot_predicted))
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='Perfect prediction')
        ax.plot([min_val, max_val], [0.9*min_val, 0.9*max_val], 'g:', alpha=0.5, label='±10%')
        ax.plot([min_val, max_val], [1.1*min_val, 1.1*max_val], 'g:', alpha=0.5)
        ax.set_xlabel('Measured ΔP (Pa)', fontsize=12)
        ax.set_ylabel('Predicted ΔP (Pa)', fontsize=12)
        ax.set_title('Hot Side (Tube)', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Calculate R² and RMSE
        hot_r2 = 1 - np.sum((np.array(hot_measured) - np.array(hot_predicted))**2) / \
                    np.sum((np.array(hot_measured) - np.mean(hot_measured))**2)
        hot_rmse = np.sqrt(np.mean((np.array(hot_measured) - np.array(hot_predicted))**2))
        hot_mape = np.mean(np.abs((np.array(hot_measured) - np.array(hot_predicted)) / np.array(hot_measured))) * 100
        
        ax.text(0.05, 0.95, f'R² = {hot_r2:.4f}\nRMSE = {hot_rmse:.1f} Pa\nMAPE = {hot_mape:.1f}%',
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"\nPlot saved to {save_path}")
        
        plt.show()
        
        return fig


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    
    optimizer = HeatExchangerOptimizer(hxs, mass_flows, pressure_changes)

    # Diagnose before optimization
    print("\n" + "="*60)
    print("BEFORE OPTIMIZATION")
    print("="*60)
    optimizer.diagnose_all()

    # Configuration matching your current tuning parameters
    param_config = {
        'K_tube_misc': (0.0, 5.0, 0.0),
        'crossflow_correction_factor': (0.5, 3.0, 1.6),
        'shell_friction_multiplier': (0.6, 2.0, 1.0),
        'nozzle_correction_factor': (0.7, 3.0, 1.0),
        'friction_divisor': (0.7, 3.0, 1.0),
        'entrance_exit_multiplier': (0.7, 2.0, 1.0),
        'K_turn': (0.0, 2.0, 1.0)
    }

    result = optimizer.optimize_parameters(param_config, weight_cold=2.0, weight_hot=1.0)

    # Diagnose after optimization
    print("\n" + "="*60)
    print("AFTER OPTIMIZATION")
    print("="*60)
    optimizer.diagnose_all()

    optimizer.plot_comparison(save_path='hx_optimization.png')