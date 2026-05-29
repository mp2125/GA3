import numpy as np
from math import exp, log
import parameters
from compressor_characteristics import k_hose_cold, k_hose_hot


class ShellAndTubeHeatExchanger:
    """
    Shell-and-tube heat exchanger model following GA3 handout correlations
    """

    def __init__(
        self,
        number_of_tubes,
        number_of_baffles,
        tube_length,
        tube_pitch=None,
        is_square_layout=False,
        tube_passes=2,
        shell_passes=1
    ):
        # Geometry
        self.number_of_tubes = number_of_tubes
        self.number_of_baffles = number_of_baffles
        self.tube_length = tube_length
        self.tube_pitch = tube_pitch
        self.is_square_layout = is_square_layout

        # constant parameters
        self.tube_outer_diameter = parameters.do
        self.tube_inner_diameter = parameters.di
        self.shell_inner_diameter = parameters.ds
        self.shell_passes = shell_passes
        self.hose_diameter = parameters.hose_diameter
        self.hose_length = parameters.hose_length
        self.tube_passes = tube_passes
        self.tube_roughness = 0.0045e-3  # approximate roughness of tube
        self.K_hose_cold = k_hose_cold  # determined from mdot max on compressor
        self.K_hose_hot = k_hose_hot

        if tube_pitch is not None:
            self.tube_pitch = tube_pitch
        else:
            phi = self.number_of_tubes * (self.tube_outer_diameter/self.shell_inner_diameter)**2
            if self.is_square_layout:
                self.tube_pitch = self.tube_outer_diameter * np.sqrt(np.pi / (4 * phi))
            else:
                self.tube_pitch = self.tube_outer_diameter * np.sqrt(np.pi / (2 * np.sqrt(3) * phi))

        # Nozzles
        self.nozzle_area_shell_side = parameters.A_noz
        self.nozzle_area_tube_side = parameters.A_noz

        # Fluid properties
        self.fluid_density = parameters.rho_w
        self.fluid_viscosity = parameters.mu

        # tuning parameters
        self.K_tube_misc = 0
        self.crossflow_correction_factor = 0.8 if not self.is_square_layout else 0.95
        # self.shell_friction_a = 0.34 if self.is_square_layout else 0.2  # Kern correlation multiplier (was 'a')
        # self.shell_friction_a *= 1
        self.hot_nozzle_correction_factor = 0.5
        self.cold_nozzle_correction_factor = 1
        self.friction_divisor = 1
        self.entrance_exit_multiplier = 1/(self.tube_passes**1.5)
        self.K_turn = 0.8 #* (self.shell_inner_diameter / self.baffle_spacing)**-2 # K ~ 2.0 for 180° turn, but reduced due to gradual turning

        self.cold_side_contributions = {
        }

    # ------------------------------------------------------------
    # GEOMETRY HELPERS
    # ------------------------------------------------------------

    @property
    def baffle_spacing(self):
        lengthCorrection = (self.number_of_baffles+2)*0.0015
        return (self.tube_length - lengthCorrection) / (self.number_of_baffles + 1)

    @property
    def crossflow_area_shell_side(self):
        """
        Minimum crossflow area between baffles (Kern-style)
        """
        return (
            self.shell_inner_diameter
            / self.shell_passes
            * self.baffle_spacing
            * (self.tube_pitch - self.tube_outer_diameter)
            / self.tube_pitch
            / self.crossflow_correction_factor
        )

    @property
    def tube_side_flow_area(self):
        """Flow area based on INNER diameter"""
        return (self.number_of_tubes / self.tube_passes) * (np.pi * self.tube_inner_diameter**2 / 4)
    
    @property
    def tubesheet_area(self):
        """Total frontal area of tubesheet"""
        return np.pi * self.shell_inner_diameter**2 / 4
    
    @property
    def sigma(self):
        """Area ratio for entrance/exit loss coefficients (Figure 8)"""
        tubes_per_pass = self.number_of_tubes / self.tube_passes
        tube_area = tubes_per_pass * np.pi * self.tube_inner_diameter**2 / 4
        return tube_area / self.tubesheet_area
    
    @property
    def shell_side_equivalent_diameter(self):

        Pt = self.tube_pitch
        do = self.tube_outer_diameter

        if self.is_square_layout:

            return (
                4 * (Pt**2 - np.pi * do**2 / 4)
                / (np.pi * do)
            )

        else:
            # triangular pitch

            return (
                4 * ((np.sqrt(3)/4) * Pt**2 - np.pi * do**2 / 8)
                / (np.pi * do / 2)
            )

    # ------------------------------------------------------------
    # ENTRANCE/EXIT LOSS COEFFICIENTS
    # ------------------------------------------------------------

    @staticmethod
    def kc_lam(L_over_D, Re, sigma):
        """Laminar entrance loss coefficient"""
        if 4 * (L_over_D) / Re < 0.05:
            print("Warning: 4L/D / Re below 0.05")
        x = 4 * (L_over_D) / Re
        offset = 1.19 * exp(-10.65 * x**0.597)
        kc = (1.08 - 0.41 * sigma) - offset
        return kc

    @staticmethod
    def kc_turb(Re, sigma):
        """Turbulent entrance loss coefficient (Re > 3000)"""
        offset = 0.14 * (1 - exp(-0.00136 * (Re - 3000)**0.622))
        kc = (0.54 - 0.39 * sigma) - offset
        return kc

    @staticmethod
    def ke_turb(Re, sigma):
        """Turbulent exit loss coefficient"""
        offset = 0.12 / (1 + 0.42 * (log(Re / 3000))**2.05)
        ke = (1 + 0.85 * sigma) * (1 - sigma)**2.35 - sigma * offset
        return ke

    @staticmethod
    def ke_lam(L_over_D, Re, sigma):
        """Laminar exit loss coefficient"""
        if 4 * (L_over_D) / Re < 0.05:
            print("Warning: 4L/D / Re below 0.05")
        x = 4 * (L_over_D) / Re
        offset = 1.19 * (exp(-10.65 * x**0.597))
        ke = (1 - sigma)**2.18 * (1 + 1.05 * sigma - 0.62 * sigma**2) - 0.66 * sigma + sigma * offset
        return ke

    def get_entrance_exit_coefficients(self, reynolds_number):
        """
        Get Kc and Ke based on flow regime and geometry
        """
        L_over_D = self.tube_length / self.tube_inner_diameter
        sigma = self.sigma
        
        # Determine flow regime
        if reynolds_number < 2300:
            # Laminar flow
            Kc = self.kc_lam(L_over_D, reynolds_number, sigma)
            Ke = self.ke_lam(L_over_D, reynolds_number, sigma)
        elif reynolds_number < 3000:
            # Transition region - interpolate
            Kc_lam = self.kc_lam(L_over_D, 2300, sigma)
            Kc_turb = self.kc_turb(3000, sigma)
            Kc = Kc_lam + (Kc_turb - Kc_lam) * (reynolds_number - 2300) / 700
            
            Ke_lam = self.ke_lam(L_over_D, 2300, sigma)
            Ke_turb = self.ke_turb(3000, sigma)
            Ke = Ke_lam + (Ke_turb - Ke_lam) * (reynolds_number - 2300) / 700
        else:
            # Turbulent flow
            Kc = self.kc_turb(reynolds_number, sigma)
            Ke = self.ke_turb(reynolds_number, sigma)
            
        return Kc*self.entrance_exit_multiplier, Ke*self.entrance_exit_multiplier
    
    def friction_factor(self, reynolds_number, relative_roughness):
        """
        Friction factor from Moody diagram (Figure 7)
        Uses Colebrook-White equation for turbulent flow
        """
        if reynolds_number < 2300:
            # Laminar
            return 64 / reynolds_number
        else:
            # Turbulent - Colebrook-White (implicit)
            # Simplified using Swamee-Jain explicit approximation
            
            # if relative_roughness < 1e-6:
            #     # Smooth tube (Blasius)
            #     return 0.316 * reynolds_number**(-0.25)
            # else:
            #     # Rough tube
            #     term1 = relative_roughness / 3.7
            #     term2 = 5.74 / (reynolds_number**0.9)
            #     return 0.25 / (np.log10(term1 + term2)**2)

            return (1.82*np.log10(reynolds_number) - 1.64)**(-2)


    # ------------------------------------------------------------
    # HOSE LOSSES (applied to BOTH inlet and outlet)
    # ------------------------------------------------------------

    def hose_velocity(self, mass_flow_rate):
        area = np.pi * self.hose_diameter**2 / 4
        return mass_flow_rate / (self.fluid_density * area)
    
    def hose_reynolds_number(self, mass_flow_rate):
        return self.hose_velocity(mass_flow_rate) * self.hose_diameter * self.fluid_density / self.fluid_viscosity
    
    def hose_pressure_drop(self, mass_flow_rate, K_hose):
        """
        Pressure drop in ONE hose using specified loss coefficient K_hose.
        Total system has 2 hoses per side (handled elsewhere).
        """
        v = self.hose_velocity(mass_flow_rate)

        return K_hose * 0.5 * self.fluid_density * v**2
    # ------------------------------------------------------------
    # SHELL SIDE (COLD FLUID)
    # ------------------------------------------------------------

    def shell_side_velocity(self, mass_flow_rate_cold):
        """Characteristic velocity through tube bundle"""
        area = self.crossflow_area_shell_side
        return (mass_flow_rate_cold / (self.fluid_density * area))

    def shell_side_reynolds_number(self, velocity):
        return (self.fluid_density * velocity * self.shell_side_equivalent_diameter) / self.fluid_viscosity

    def shell_side_row_drag_coefficient(self, Re):

        pitch_ratio = (
            self.tube_pitch
            / self.tube_outer_diameter
        )

        if self.is_square_layout:

            return (
                2.1
                * Re**(-0.2)
                * pitch_ratio**(-0.5)
            )

        else:

            return (
                # 1.8
                # * Re**(-0.15)
                # * pitch_ratio**(-0.4)
                1.8
                * Re**(-0.15)
                * pitch_ratio**(-0.4)
            )

    def shell_side_pressure_drop(self, mass_flow_rate_cold):
        """
        Total shell-side pressure drop with:
        - Bundle crossflow (eq 9): ΔP = 4 * a * Re^(-0.15) * N * rho * V^2
        - Window zone turning losses at each baffle
        - Nozzle losses: 2 dynamic heads
        - Multiplied by shell_passes for multi-pass configurations
        """
        velocity = self.shell_side_velocity(mass_flow_rate_cold)
        reynolds = self.shell_side_reynolds_number(velocity)
        # a = self.shell_side_friction_coefficient(reynolds)
        
        # # Number of tube rows crossed per shell pass (approximation)
        # N = self.number_of_baffles + 1
        
        # # Bundle pressure drop (equation 9) - per shell pass
        # bundle_pressure_drop_per_pass = (
        #     4 * a * N * self.fluid_density * velocity**2
        # )

        # f = self.shell_side_friction_coefficient(reynolds)
        rows_per_section = self.shell_inner_diameter / self.tube_pitch
        total_rows = rows_per_section * (self.number_of_baffles+1)

        zeta = self.shell_side_row_drag_coefficient(reynolds)

        bundle_pressure_drop_per_pass = (
            zeta
            * rows_per_section
            * total_rows
            * 0.5
            * self.fluid_density
            * velocity**2
        )

        total_bundle_pressure_drop = bundle_pressure_drop_per_pass*self.shell_passes
        self.cold_side_contributions['bundle'] = total_bundle_pressure_drop

        
        # Window zone turning losses - 180° turn at each baffle
        turning_losses_per_pass = (
            self.number_of_baffles  # one turn per baffle
            * self.K_turn
            * 0.5
            * self.fluid_density
            * velocity**2
        )

        total_turning_losses = turning_losses_per_pass * self.shell_passes
        end_turning_loss = 2 * 0.5 * self.fluid_density * velocity**2 * (self.shell_passes-1)


        self.cold_side_contributions['turning'] = total_turning_losses
        
        # Total bundle pressure drop accounting for shell passes
        bundle_pressure_drop = total_bundle_pressure_drop + total_turning_losses + end_turning_loss
        
        # Nozzle losses: 2 dynamic heads (inlet and outlet only, not per pass)
        nozzle_velocity = mass_flow_rate_cold / (
            self.fluid_density * self.nozzle_area_shell_side
        ) * 1
        nozzle_pressure_drop = 2 * self.cold_nozzle_correction_factor * 0.5 * self.fluid_density * nozzle_velocity**2

        self.cold_side_contributions['nozzle'] = nozzle_pressure_drop
        
        return  bundle_pressure_drop + nozzle_pressure_drop

    # ------------------------------------------------------------
    # TUBE SIDE (HOT FLUID)
    # ------------------------------------------------------------

    def tube_side_velocity(self, mass_flow_rate_hot):
        """Velocity based on INNER diameter"""
        area = self.tube_side_flow_area
        return mass_flow_rate_hot / (self.fluid_density * area)

    def tube_side_reynolds_number(self, velocity):
        """Reynolds number based on INNER diameter"""
        return (self.fluid_density * velocity * self.tube_inner_diameter) / self.fluid_viscosity


    def tube_side_pressure_drop(self, mass_flow_rate_hot):
        """
        Total tube-side pressure drop:
        1. Friction losses in tubes (Moody diagram)
        2. Entrance/exit losses (equation 8, Figure 8)
        3. Nozzle losses (2 dynamic heads)
        4. Misc Losses
        """
        velocity = self.tube_side_velocity(mass_flow_rate_hot)
        reynolds = self.tube_side_reynolds_number(velocity)
        relative_roughness = self.tube_roughness / self.tube_inner_diameter
        friction_factor = self.friction_factor(reynolds, relative_roughness) / self.friction_divisor

        # 1. Pipe friction (Darcy-Weisbach): ΔP = f * (L/D) * ρ * V^2 / 2
        effective_length = self.tube_length * self.tube_passes
        pipe_friction_loss = (
            friction_factor
            * (effective_length / self.tube_inner_diameter)
            * 0.5
            * self.fluid_density
            * velocity**2
        )

        self.pipe_friction_loss = pipe_friction_loss

        # 2. Entrance/exit losses (equation 8): ΔP = 0.5 * ρ * V^2 * (Kc + Ke)
        Kc, Ke = self.get_entrance_exit_coefficients(reynolds)
        # Multiply by tube_passes since entrance/exit occurs at each pass
        entrance_exit_loss = (
            self.tube_passes
            * 0.5
            * self.fluid_density
            * velocity**2
            * (Kc + Ke)
        )

        self.entrance_exit_loss = entrance_exit_loss

        # 3. Nozzle losses: 2 dynamic heads
        nozzle_velocity = mass_flow_rate_hot / (
            self.fluid_density * self.nozzle_area_tube_side
        )
        self.nozzle_velocity_hot = nozzle_velocity
        nozzle_loss = 2 * self.hot_nozzle_correction_factor * 0.5 * self.fluid_density * nozzle_velocity**2

        turning_loss = 2 * 0.5 * self.fluid_density * velocity**2 * (self.tube_passes-1)


        # 4. Misc Losses
        misc_loss = (
            self.K_tube_misc
            * 0.5
            * self.fluid_density
            * velocity**2
        )

        # pipe_friction_loss = 0
        return pipe_friction_loss + entrance_exit_loss + nozzle_loss + turning_loss + misc_loss

    # ------------------------------------------------------------
    # SYSTEM INTERFACE
    # ------------------------------------------------------------

    def cold_side_pressure_drop(self, mass_flow_rate_cold):
        exchanger_dp = self.shell_side_pressure_drop(mass_flow_rate_cold)

        hose_dp = 2 * self.hose_pressure_drop(
            mass_flow_rate_cold,
            self.K_hose_cold
        )

        self.cold_side_contributions['hose'] = hose_dp

        return exchanger_dp + hose_dp

    def hot_side_pressure_drop(self, mass_flow_rate_hot):
        exchanger_dp = self.tube_side_pressure_drop(mass_flow_rate_hot)

        hose_dp = 2 * self.hose_pressure_drop(
            mass_flow_rate_hot,
            self.K_hose_hot
        )

        return exchanger_dp + hose_dp

    def total_pressure_drop(self, mass_flow_rate_cold, mass_flow_rate_hot):
        return (
            self.cold_side_pressure_drop(mass_flow_rate_cold),
            self.hot_side_pressure_drop(mass_flow_rate_hot),
        )